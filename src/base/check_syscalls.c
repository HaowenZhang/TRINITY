#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <inttypes.h>
#include <sys/socket.h>
#include <sys/wait.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <sys/errno.h>
#include <unistd.h>
#include <signal.h>
#include "inthash.h"
#include "check_syscalls.h"

//#define DEBUG_IO

char *unread = NULL;
int64_t unread_size = 0;

void system_error(char *errmsg) {
  fprintf(stderr, "[Error] %s\n", errmsg);
  perror("[Error] Reason");
  exit(1);
}

pid_t check_waitpid(pid_t pid) {
  int stat_loc;
  pid_t res;
  do {
    res = waitpid(pid, &stat_loc, 0);
  } while ((res < 0) && (errno == EINTR));

  if (res < 0) system_error("Waiting for child process failed.");
  return res;
}


FILE *check_fopen(char *filename, char *mode) {
  FILE *res = fopen(filename, mode);
  if (res == NULL) {
    if (mode[0] == 'w')
      fprintf(stderr, "[Error] Failed to open file %s for writing!\n", filename);
    else if (mode[0] == 'a')
      fprintf(stderr, "[Error] Failed to open file %s for appending!\n", filename);
    else
      fprintf(stderr, "[Error] Failed to open file %s for reading!\n", filename);
    exit(1);
  }
#ifdef DEBUG_IO
  fprintf(stderr, "[Note] Opened file %s with mode '%s' in fileno %d.\n", 
	  filename, mode, fileno(res));
#endif /* DEBUG_IO */
  return res;
}

FILE *check_popen(char *command, char *mode) {
  FILE *res = popen(command, mode);
  if (res == NULL) {
    fprintf(stderr, "[Error] Failed to start command %s!\n", command);
    exit(1);
  }
#ifdef DEBUG_IO
  fprintf(stderr, "[Note] Opened command %s with mode '%s' in fileno %d.\n", 
	  command, mode, fileno(res));
#endif /* DEBUG_IO */
  return res;
}

FILE *check_rw_socket(char *command, pid_t *pid) {
  int sockets[2], status;
  pid_t wres;
  FILE *res;
  if (socketpair(AF_UNIX, SOCK_STREAM, 0, sockets)<0)
    system_error("Failed to create socket pair!");
  *pid = fork();
  if (*pid < 0) system_error("Failed to fork new process!");
  if (!*pid) {
    if (dup2(sockets[1], 0) < 0) system_error("Failed to reopen stdin!");
    if (dup2(sockets[1], 1) < 0) system_error("Failed to reopen stdout!");
    close(sockets[0]);
    if (execlp("sh", "sh", "-c", command, NULL) < 0)
      system_error("Failed to exec command!");
  }
  close(sockets[1]);
  res = fdopen(sockets[0], "r+");
  if (!res) system_error("Failed to convert socket to stream!");
  do {
    wres = waitpid(*pid, &status, WNOHANG);
  } while ((wres < 0) && (errno == EINTR));
  if (wres < 0) {
    fprintf(stderr, "[Error] Failed to start child process: %s\n", command);
    exit(1);
  }
#ifdef DEBUG_IO
  fprintf(stderr, "[Note] Started command %s with mode 'r+' in fileno %d.\n", 
	  command, fileno(res));
#endif /* DEBUG_IO */
  return res;
}

void check_lseek(int fd, off_t offset, int whence) {
  int64_t res = lseek(fd, offset, whence);
  if (res<0) {
    fprintf(stderr, "[Error] Lseek error in fileno %d: ", fd);
    perror(NULL);
    exit(1);
  }
}

void rw_socket_close(FILE *res, pid_t pid) {
  fclose(res);
  kill(pid, 9);
  check_waitpid(pid);
}

void *check_realloc(void *ptr, size_t size, char *reason) {
  if (size > 0) {
    void *res = realloc(ptr, size);
    if (res == NULL) {
      fprintf(stderr, "[Error] Failed to allocate memory (%s)!\n", reason);
      exit(1);
    }
    return res;
  }
  if (ptr != NULL) free(ptr);
  return(NULL);
}


void _io_err(int rw, size_t size, size_t nitems, FILE *stream) {
  char *verb = (rw) ? "write" : "read";
  char *dir = (rw) ? "to" : "from";
  char *items = (nitems == 1) ? "item" : "items";

  fprintf(stderr, "[Error] Failed to %s %"PRIu64" %s of size "
	  "%"PRIu64" bytes %s fileno %d!\n", 
	  verb, (uint64_t)nitems, items, (uint64_t)size, dir, fileno(stream));
  if (feof(stream))
    fprintf(stderr, "[Error] Reason: end of file (offset %"PRIu64").\n",
	    (uint64_t)ftello(stream));
  else
    perror("[Error] Reason");
  exit(1);
}

void check_fseeko(FILE *stream, off_t offset, int whence) {
  if (fseeko(stream, offset, whence) < 0) {
    fprintf(stderr, "[Error] Seek error in fileno %d: ", fileno(stream));
    perror(NULL);
    exit(1);
  }
}

//Works even for pipes
void check_fskip(FILE *stream, off_t offset, char *buffer, size_t buf_size) {
  int64_t n = 0;
  while (n<offset) {
    int64_t to_read = offset-n;
    if (buf_size < to_read) to_read = buf_size;
    n += check_fread(buffer, 1, to_read, stream);
  }
}

void check_limited_funread(void *ptr, size_t size, size_t nitems) {
  if (unread_size) {
    fprintf(stderr, "[Error] Tried to unread twice in a row\n");
    exit(1);
  }
  check_realloc_s(unread, size, nitems);
  unread_size = size*nitems;
  memcpy(unread, ptr, unread_size);  
}

size_t check_fread(void *ptr, size_t size, size_t nitems, FILE *stream) {
  size_t res = 1, nread = 0;
  if (unread_size) {
    if (unread_size != (size*nitems)) {
      fprintf(stderr, "[Error] funread must be followed by identical fread!\n");
      exit(1);
    }
    memcpy(ptr, unread, unread_size);
    check_realloc_s(unread, 0, 0);
    unread_size = 0;
    return nitems;
  }

  while (nread < nitems) {
    res = fread(ptr, size, nitems-nread, stream);
    if (res <= 0) _io_err(0, size, nitems, stream);
    nread += res;
    ptr = ((char *)ptr) + res*size;
  }
  return nread;
}

char *check_fgets(char *ptr, size_t size, FILE *stream) {
  char *res = fgets(ptr, size, stream);
  if (!res) _io_err(0, size, 1, stream);
  return res;
}


size_t check_fwrite(void *ptr, size_t size, size_t nitems, FILE *stream) {
  size_t res = 1, nwritten = 0;
  while (nwritten < nitems) {
    res = fwrite(ptr, size, nitems, stream);
    if (res <= 0) _io_err(1, size-1, nitems, stream);
    nwritten += res;
  }
  return nwritten;
}

void *check_mmap_file(char *filename, char mode, int64_t *length) {
  FILE *tf;
  int flags = MAP_SHARED, prot = PROT_READ;
  struct stat ts;
  if (mode == 'r') tf = check_fopen(filename, "rb");
  else if (mode == 'w') {
    tf = check_fopen(filename, "r+b");
    prot |= PROT_WRITE;
  }
  else {
    fprintf(stderr, "[Error] Invalid mode %c passed to check_mmap_file!\n", mode);
    exit(1);
  }
  int fd = fileno(tf);
  if (fstat(fd, &ts)!=0) {
    fprintf(stderr, "[Error] Fstat failure on file %s!\n", filename);
    perror("[Error] Reason");
    exit(1);
  }
  void *res = mmap(NULL, ts.st_size, prot, flags, fd, 0);
  if (res == MAP_FAILED) {
    fprintf(stderr, "[Error] Mmap failure on file %s, mode %c!\n", filename, mode);
    perror("[Error] Reason");
    exit(1);
  }
  fclose(tf);
  if (length) *length = ts.st_size;
  return res;
}
