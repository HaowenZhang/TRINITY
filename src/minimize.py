import os, sys
import subprocess as sp
import numpy as np
from scipy.optimize import minimize
# import python.ReadConfig as rc

# # Applies to UniverseMachine but not to Trinity.
# if len(sys.argv) < 2:
#     print("Usage:",sys.argv[0],"um.cfg")
#     sys.exit()

np.set_printoptions(linewidth=100000)

cmd = "./all_smf_mcmc 1 1 1 1 1 0 0 ../mf_bolshoi_planck.dat ../obs_AGN_old_ABHMF/*.* ../obs_AGN_old_ABHMF/Ueda_QLF_CTK/*.qlf ../obs_AGN_old_ABHMF/Aird_qpdf/*no_highest_0.3dex.qpdf ../obs_AGN_old_ABHMF/qf/*.qf ../obs_AGN_old_ABHMF/uvlf/new/*.uvlf"

p = sp.Popen(cmd, shell=True, bufsize=1024,
             stdin=sp.PIPE, stdout=sp.PIPE, close_fds=True)

#steps = np.array([0.1, 0.2, 0.2, 0.2, 0.1, 0.1, 0.1, 0.3, 0.3, 1.0, 1.0, 1.0, 1.0, 0.02, 0.1, 0.2, 0.2, 1, 1, 1, 0.1, 0.02, 0.05, 0.02, 0.03, 0.01, 0.02, 0.003, 0.002, 1, 1, 1, 1, 1])

# eff_0, 1, 2, 3
# v_0, 1, 2, 3
# alpha_0, 1, 2, 3
# beta_0, 1, 2
# mu
# kappa
# scatter
# sigma_z
# icl_frac
# mu_a
# scatter_a
# rho_05

# bh_beta_0, 1, 2
# bh_gamma_0, 1, 2
# bh_alpha_0, 1
# bh_delta_0
# bh_duty_0, 1
# bh_merge_f_0, 1
# bh_scatter_1
# dc_mbh_w_0, 1
# eta_mu_0
# eff_0
# bh_delta_0

steps = np.array([0.01, 0.03, 0.03, 0.03, \
                  0.1, 0.3, 0.3, 0.3,\
                  0.1, 0.3, 0.3, 0.3,\
                  0.1, 0.3, 0.3,\
                  0.01,\
                  0.03,\
                  0.03,\
                  0.01,\
                  0.1,\
                  0.01,\
                  0.01,\
                  0.01,\

                  0.1,0.03,0.01,\
                  0.1,0.03,0.01,\
                  0.1, 0.1,\
                  0.1,\
                  0.1, 0.1,\
                  0.1, 0.1,\

                  0.03,\
                  0.3, 0.1,\
                  0.1,\
                  0.03,\
                  0.1,\

                  0.3, 0.3, 0.03,\
                  0.1, 0.1, 0.01])

chi2_dict = {}

def get_chi2(x):
#    global chi2_dict
#    dkey = str(x)
#    if (dkey in chi2_dict):
#        return chi2_dict[dkey]
    x = np.multiply(x, steps)
    # x = np.concatenate((x[0:26], np.array([1]), x[26:36], np.array([0.07]), x[36:45]))
    x = np.concatenate((x[0:4], np.array([4]), x[4:15], np.array([0.055]),\
        np.zeros(8), x[15:21], np.array([0]), np.atleast_1d(x[21]), np.array([0]), np.array([0,0,1]), x[22:36],\
        np.array([0,0,-1.5,0,0.3]), np.atleast_1d(x[36]), np.array([0.2]), x[37:48]))

    y = str(str(x)[1:-1])
    print(y)
    p.stdin.write(str.encode(y+" 0\n"))
    p.stdin.flush()
    chi2 = float(p.stdout.readline())
#    if (chi2>1e20):
#        chi2 = -1.
    print(y," ",chi2)
    sys.stdout.flush()
#    chi2_dict[dkey] = chi2
    return chi2


# x = np.array([0.111566383751, 1.851165306043, -0.262498814697, 2.453342000000, 2.060578524077, 0.322385485956, -0.095760500000, 0.706293003179, -6.678835123879, -0.742246653101, -0.268311073596, 2.412981052222, -1.103871743392, 0.382719417322, 0.871704118926, 0.080322851043, -1.313790240545, 3.342866290000, -0.972355917846, -0.169145118576, 1.674300061812, 0.290904124799, 2.827757777301, -0.569456960882, 1.896251348213, -0.080875527084, 2.322659214735, -0.144645720025, 0.210928240000, 0.466701980000, -1.126892182934, 1.080337094572, 0.081656800292, -0.125073476035, 0.322800227535, 0.059910483272, -0.170076796082, -0.108142107758, 0.315259793596, 9.103332360805, -0.040651324704, 0.442313450000, 0.560758442785, 0.472760328977])
#x = np.array([0.103966668408056,1.83063455621321,-0.277654271380385,2.48415368092332,2.06070261079063,0.321276981182992,-0.099259307971709,0.698611015909492,-6.78203027787717,-0.739069324146649,-0.242721698327971,2.35968078709715,-0.769998322065594,-0.116483508267006,0.714846466734106,0.0782099967020499,-1.30935957816679,3.34888593019753,-0.968024308196845,-0.350872870921099,1.62602563488957,0.422315986805749,3.05486114746315,0.362931710371528,2.95574055440488,-0.109709699148675,2.32098601411566,-0.155812123549582,0.208789379382599,0.451639335282312,-1.13797216523246,1.06671891715612,0.102255261540046,-0.230386628725287,0.302581799836008,0.0393551448374071,-0.16244070603523,-0.125857572293462,0.330001386022831,8.97192979696666,-0.0440697344531944,0.431915895223722,0.600482034644315, 0.600482034644315])
x = np.array([0.416521, 0.795884, 1.949462, -0.010836,\
             2.237558, 1.493246, 1.320093, -0.121708,\
            -4.390630, 29.570738, 18.612680, -1.925817,\
             1.203754, 5.990015, 2.019021,\
             -0.008821, 0.270311, 0.277690,\
              0.054616, -0.431619, 0.157875,\
              -0.027740,\
              0.615966,\
               8.610493, 0.087433, -0.160713,\
                1.034246, 0.443126, 0.008288,\
                -0.277316, 1.044703, 2.013441,\
                 0.888630, -0.029875,\
                  -0.014546, 0.030973,\
                  -0.342050,\
                  8.396402, 0.583703,\
                   0.135861,\
                    -1.281819,\
                     1.219345,\
                      12.280687, 0.817332, 0.635136,\
                       0.562134, 0.673732, 0.152052])

print(len(x), len(steps))

# x = np.concatenate((x[0:26], np.array([1]), x[26:36], np.array([0.07]), x[36:45]))

# print(x[0:4].shape)
# print(np.array([4]).shape)
# print(x[4:15].shape)
# print(np.array([0.055]).shape)

# print(np.zeros(8).shape)
# print(x[15:21].shape)
# print(np.array([0]).shape)
# print(np.array(x[21]).shape)

# x = np.concatenate((x[0:4], np.array([4]), x[4:15], np.array([0.055]),\
#         np.zeros(8), x[15:21], np.array([0]), np.atleast_1d(x[21]), np.array([0]), np.array([0,0,1]), x[22:36],\
#         np.array([0,0,-1.5,0,0.3]), np.atleast_1d(x[36]), np.array([0.2]), x[37:48]))





print(len(x))
print(x)

x = np.divide(x, steps)
res = minimize(get_chi2, x, method='COBYLA', options={'rhobeg':1e-1, 'maxiter': 10000}) #, constraints=({'type':'ineq', 'fun':get_chi2}))
chi2 = get_chi2(res.x)
x = np.multiply(res.x, steps)
x = np.concatenate((x[0:26], np.array([1]), x[26:36], np.array([0.07]), x[36:45]))
print("Final answer:")
print(str(str(x)[1:-1])," ",chi2)

p.terminate()

# cfg = rc.ReadConfig(sys.argv[1])
# if "OUTBASE" in cfg.data:
#     file = open(cfg.data["OUTBASE"]+"/"+"done", "w")
#     file.close()


