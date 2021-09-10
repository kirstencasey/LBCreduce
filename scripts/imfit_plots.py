import lbcred, glob, os, yaml
import matplotlib.pyplot as plt
import numpy as np
from lbcred.utils import io
from lbcred.model import imfit

# Define files we want to read
in_dir = '/Users/kirstencasey/real_image_tests/imfit_sbf/'
save_dir = '/Users/kirstencasey/real_image_tests/imfit_sbf/'
config_filename = os.path.join('../modeling-config_imfit.yml')
fn_stub = 'real_image_tests_fixed_artpop_fix_sersic_on_b_autogenmask'
fn_id_r = f'cutout_SEsky_mock_injected_Sersic_bestfit-params_{fn_stub}_iter*_r.txt'
fn_id_b = f'cutout_SEsky_mock_injected_Sersic_bestfit-params_{fn_stub}_iter*_b.txt'
files_r = glob.glob(os.path.join(in_dir,fn_id_r))
files_b = glob.glob(os.path.join(in_dir,fn_id_b))
random_artpops = True
random_sersics = False
model_funcs = ['Sersic']

# Open and read config file
with open(config_filename, 'r') as filename:
    config = yaml.load(filename, Loader=yaml.FullLoader)

def sort_key(fn):
    return int(fn.split('_iter')[-1].split('_')[0])
files_r.sort(key=sort_key)
files_b.sort(key=sort_key)

# Read files
chisq_r = []
chisq_b = []
radii = []
pa = []
ell = []
n = []
I_e_r = []
mags_r = []
xpos = []
ypos = []
I_e_b = []
mags_b = []
iters = []
sky_r = []
sky_b = []
x_slopes_r = []
y_slopes_r = []
x_slopes_b = []
y_slopes_b = []
num_iters = len(files_r)
comp = 'comp_1'
r_true = 149.323 #221.02 #224.133
pa_true = 0 #9.53788
ell_true = 0 #-0.0926788
n_true = 0.52078
I_e_r_true = 89.7946 # 8.4 # 89.7946
mag_r_true = 15.767578 #14.76083716950836 # 14.531013184162031 # 14.594793336161725
xpos_true = 525.0
ypos_true = 525.0
I_e_b_true = 23.825 # 2.825 # 23.825
mag_b_true = 17.030594 #16.019337153610394 # 15.780570355571324 # 15.861027188510267
color_true =  mag_b_true - mag_r_true
sbf_mag_true = 26.70403490761306 # 26.734176700344815
true_dist = 3.7 # Mpc

plot_I_e = False
if random_artpops: model_type = 'artpop'
elif random_sersics : model_type = 'sersic'

if random_artpops or random_sersics:
    plot_I_e = False
    if not random_sersics:
        file = open(in_dir+f'{model_type}_parameters_{fn_stub}_num-iters{num_iters}_true_sbf_magnitude', "rb")
        sbf_mag_true = np.load(file,allow_pickle=True)
        file.close

    file = open(in_dir+f'{model_type}_parameters_{fn_stub}_num-iters{num_iters}_n', "rb")
    n_true = np.load(file,allow_pickle=True)
    file.close

    file = open(in_dir+f'{model_type}_parameters_{fn_stub}_num-iters{num_iters}_theta', "rb")
    pa_true = np.load(file,allow_pickle=True)
    file.close

    file = open(in_dir+f'{model_type}_parameters_{fn_stub}_num-iters{num_iters}_ellipticity', "rb")
    ell_true = np.load(file,allow_pickle=True)
    file.close

    file = open(in_dir+f'{model_type}_parameters_{fn_stub}_num-iters{num_iters}_distance', "rb")
    true_dist = np.load(file,allow_pickle=True)
    file.close

    file = open(in_dir+f'{model_type}_parameters_{fn_stub}_num-iters{num_iters}_radius', "rb")
    r_true = np.load(file,allow_pickle=True)
    if random_artpops:
        r_true = np.arctan(r_true/(true_dist*1e6)) * 180 / np.pi * 3600 / config['pixscale'] # convert to pixels
    file.close

    file = open(in_dir+f'{model_type}_parameters_{fn_stub}_num-iters{num_iters}_mag_r', "rb")
    mag_r_true = np.load(file,allow_pickle=True)
    file.close

    file = open(in_dir+f'{model_type}_parameters_{fn_stub}_num-iters{num_iters}_mag_b', "rb")
    mag_b_true = np.load(file,allow_pickle=True)
    file.close

    color_true =  mag_b_true - mag_r_true

for i in range(len(files_r)):
    fn_r = files_r[i]
    fn_b = files_b[i]

    bestfit_r = io.read_results(fn_r, model_funcs )
    bestfit_b = io.read_results(fn_b, model_funcs )

    chisq_r.append(bestfit_r['reduced_chisq'])
    chisq_b.append(bestfit_b['reduced_chisq'])
    radii.append(bestfit_r[comp]['r_e'])
    pa.append(bestfit_r[comp]['PA'])
    ell.append(bestfit_r[comp]['ell'])
    n.append(bestfit_r[comp]['n'])
    xpos.append(bestfit_r[comp]['X0'])
    ypos.append(bestfit_r[comp]['Y0'])
    I_e_r.append(bestfit_r[comp]['I_e'])
    I_e_b.append(bestfit_b[comp]['I_e'])
    if comp == 'comp_2' and model_funcs[0]=='TiltedSkyPlane':
        sky_r.append(bestfit_r['comp_1']['I_0'])
        x_slopes_r.append(bestfit_r['comp_1']['m_x'])
        y_slopes_r.append(bestfit_r['comp_1']['m_y'])
        sky_b.append(bestfit_b['comp_1']['I_0'])
        x_slopes_b.append(bestfit_b['comp_1']['m_x'])
        y_slopes_b.append(bestfit_b['comp_1']['m_y'])
    elif comp == 'comp_2' and model_funcs[0]=='FlatSky':
        sky_r.append(bestfit_r['comp_1']['I_sky'])
        sky_b.append(bestfit_b['comp_1']['I_sky'])

    iters.append(i)

    mag_r, mag_b, color = imfit.summarize_results(config, bestfit_r[comp], bestfit_b[comp])
    mags_r.append(mag_r)
    mags_b.append(mag_b)

chisq_r=np.asarray(chisq_r)
chisq_b=np.asarray(chisq_b)
radii=np.asarray(radii)
pa=np.asarray(pa)
ell=np.asarray(ell)
n=np.asarray(n)
xpos=np.asarray(xpos)
ypos=np.asarray(ypos)
I_e_r=np.asarray(I_e_r)
I_e_b=np.asarray(I_e_b)
mags_r=np.asarray(mags_r)
mags_b=np.asarray(mags_b)
iters=np.asarray(iters)
sky_r=np.asarray(sky_r)
sky_b=np.asarray(sky_b)
x_slopes_r=np.asarray(x_slopes_r)
y_slopes_r=np.asarray(y_slopes_r)
x_slopes_b=np.asarray(x_slopes_b)
y_slopes_b=np.asarray(y_slopes_b)
colors=mags_b-mags_r
print('r_true:', r_true, '\nradii: ', radii)
radii_err = (radii - r_true) / r_true
pa_err = pa - pa_true
ell_err = ell - ell_true
n_err = (n - n_true) / n_true
xpos_err = xpos - 1 - xpos_true
ypos_err = ypos - 1 - ypos_true
I_e_r_err = (I_e_r - I_e_r_true) / I_e_r_true
I_e_b_err =(I_e_b - I_e_b_true) / I_e_b_true
mags_r_err = mags_r - mag_r_true
mags_b_err = mags_b - mag_b_true
colors_err = colors - color_true
pos_err = np.sqrt(xpos_err**2 + ypos_err**2)

plt.clf()
plt.scatter(radii_err, n_err,color='purple')
plt.xlabel('Fractional Radius Error')
plt.ylabel('Sersic index Error (b-band)')
plt.axvline(x=0, color='black', linestyle='dashed', linewidth=3)
plt.axhline(y=0, color='black', linestyle='dashed', linewidth=3)
plt.savefig(save_dir+f'imfit_results_{fn_stub}_compare_radius_err_n_err.png')

if plot_I_e:
    plt.clf()
    plt.scatter(radii_err, I_e_r_err,color='purple')
    plt.xlabel('Fractional Radius Error')
    plt.ylabel('Ie Error (r-band)')
    plt.axvline(x=0, color='black', linestyle='dashed', linewidth=3)
    plt.axhline(y=0, color='black', linestyle='dashed', linewidth=3)
    plt.savefig(save_dir+f'imfit_results_{fn_stub}_compare_radius_err_I_e_r_err.png')

    plt.clf()
    plt.scatter(radii_err, I_e_b_err,color='purple')
    plt.xlabel('Fractional Radius Error')
    plt.ylabel('Ie Error (b-band)')
    plt.axvline(x=0, color='black', linestyle='dashed', linewidth=3)
    plt.axhline(y=0, color='black', linestyle='dashed', linewidth=3)
    plt.savefig(save_dir+f'imfit_results_{fn_stub}_compare_radius_err_I_e_b_err.png')

plt.clf()
plt.scatter(radii_err, mags_r_err,color='purple')
plt.xlabel('Fractional Radius Error')
plt.ylabel('Magnitude Error (r-band)')
plt.axvline(x=0, color='black', linestyle='dashed', linewidth=3)
plt.axhline(y=0, color='black', linestyle='dashed', linewidth=3)
plt.savefig(save_dir+f'imfit_results_{fn_stub}_compare_radius_err_mag_r_err.png')

plt.clf()
plt.scatter(radii_err, mags_b_err ,color='purple')
plt.xlabel('Fractional Radius Error')
plt.ylabel('Magnitude Error (b-band)')
plt.axvline(x=0, color='black', linestyle='dashed', linewidth=3)
plt.axhline(y=0, color='black', linestyle='dashed', linewidth=3)
plt.savefig(save_dir+f'imfit_results_{fn_stub}_compare_radius_err_mag_b_err.png')

plt.clf()
plt.scatter(radii_err, colors_err,color='purple')
plt.xlabel('Fractional Radius Error')
plt.ylabel('Color Error')
plt.axvline(x=0, color='black', linestyle='dashed', linewidth=3)
plt.axhline(y=0, color='black', linestyle='dashed', linewidth=3)
plt.savefig(save_dir+f'imfit_results_{fn_stub}_compare_radius_err_color_err.png')

plt.clf()
plt.scatter(xpos_err, ypos_err,color='purple')
plt.xlabel('x-pos Error')
plt.ylabel('y-pos Error')
plt.axvline(x=0, color='black', linestyle='dashed', linewidth=3)
plt.axhline(y=0, color='black', linestyle='dashed', linewidth=3)
plt.axvline(x=0, color='black', linestyle='dashed', linewidth=3)
plt.axhline(y=0, color='black', linestyle='dashed', linewidth=3)
plt.savefig(save_dir+f'imfit_results_{fn_stub}_position_err.png')

plt.clf()
plt.scatter(pos_err,radii_err,color='purple')
plt.ylabel('Fractional Radius Error')
plt.xlabel('Position Error (pixels)')
plt.axvline(x=0, color='black', linestyle='dashed', linewidth=3)
plt.axhline(y=0, color='black', linestyle='dashed', linewidth=3)
plt.savefig(save_dir+f'imfit_results_{fn_stub}_compare_radius_err_pos_err.png')

plt.clf()
plt.scatter(radii, pos_err,color='purple')
plt.xlabel('ArtPop Radius')
plt.ylabel('Position Error (pixels)')
plt.axhline(y=0, color='black', linestyle='dashed', linewidth=3)
plt.savefig(save_dir+f'imfit_results_{fn_stub}_compare_radius_pos_err.png')

plt.clf()
plt.scatter(radii, radii_err,color='purple')
plt.xlabel('ArtPop Radius')
plt.ylabel('Fractional Radius Error')
plt.axhline(y=0, color='black', linestyle='dashed', linewidth=3)
plt.savefig(save_dir+f'imfit_results_{fn_stub}_compare_radius_fracradius_err.png')
plt.clf()
plt.scatter(radii_err,chisq_b,color='purple')
plt.ylabel('Reduced Chi Squared (b-band)')
plt.xlabel('Fractional Radius Error')
plt.savefig(save_dir+f'imfit_results_{fn_stub}_compare_chisqb_radius_err.png')
plt.clf()
plt.scatter(chisq_b, chisq_r,color='purple')
plt.xlabel('Reduced Chi Squared (b-band)')
plt.ylabel('Reduced Chi Squared (r-band)')
plt.savefig(save_dir+f'imfit_results_{fn_stub}_compare_chisqb_chisqr.png')
plt.clf()
plt.scatter(pos_err, chisq_b,color='purple')
plt.ylabel('Reduced Chi Squared (b-band)')
plt.xlabel('Position Error (pixels)')
plt.savefig(save_dir+f'imfit_results_{fn_stub}_compare_chisqb_pos_err.png')
plt.clf()
plt.scatter(colors_err,chisq_b,color='purple')
plt.xlabel('Color Error')
plt.ylabel('Reduced Chi Squared (b-band)')
plt.savefig(save_dir+f'imfit_results_{fn_stub}_compare_chisqb_color_err.png')
plt.clf()
plt.scatter(mags_r_err, chisq_r, color='purple')
plt.xlabel('Magnitude Error (r-band)')
plt.ylabel('Reduced Chi Squared (r-band)')
plt.savefig(save_dir+f'imfit_results_{fn_stub}_compare_chisqr_mag_r_err.png')
plt.clf()
plt.scatter(mags_b_err,chisq_b,color='purple')
plt.xlabel('Magnitude Error (b-band)')
plt.ylabel('Reduced Chi Squared (b-band)')
plt.savefig(save_dir+f'imfit_results_{fn_stub}_compare_chisqb_mag_b_err.png')
plt.clf()
plt.scatter(n_err, chisq_b,color='purple')
plt.xlabel('Sersic Index Error')
plt.ylabel('Reduced Chi Squared (b-band)')
plt.savefig(save_dir+f'imfit_results_{fn_stub}_compare_chisqb_n_err.png')
plt.clf()
plt.scatter(ell_err, chisq_b,color='purple')
plt.xlabel('Ellipticity Error')
plt.ylabel('Reduced Chi Squared (b-band)')
plt.savefig(save_dir+f'imfit_results_{fn_stub}_compare_chisqb_ell_err.png')
plt.clf()
plt.scatter(pa_err, chisq_b,color='purple')
plt.xlabel('Position Angle Error')
plt.ylabel('Reduced Chi Squared (b-band)')
plt.savefig(save_dir+f'imfit_results_{fn_stub}_compare_chisqb_pa_err.png')
plt.clf()
plt.scatter(radii_err,chisq_r,color='purple')
plt.ylabel('Reduced Chi Squared (r-band)')
plt.xlabel('Fractional Radius Error')
plt.savefig(save_dir+f'imfit_results_{fn_stub}_compare_chisqr_radius_err.png')
plt.clf()
plt.scatter(pos_err, chisq_r,color='purple')
plt.ylabel('Reduced Chi Squared (r-band)')
plt.xlabel('Position Error (pixels)')
plt.savefig(save_dir+f'imfit_results_{fn_stub}_compare_chisqr_pos_err.png')
plt.clf()
plt.scatter(colors_err,chisq_r,color='purple')
plt.xlabel('Color Error')
plt.ylabel('Reduced Chi Squared (r-band)')
plt.savefig(save_dir+f'imfit_results_{fn_stub}_compare_chisqr_color_err.png')
plt.clf()
plt.scatter(mags_r_err, chisq_b, color='purple')
plt.xlabel('Magnitude Error (b-band)')
plt.ylabel('Reduced Chi Squared (b-band)')
plt.savefig(save_dir+f'imfit_results_{fn_stub}_compare_chisqb_mag_r_err.png')
plt.clf()
plt.scatter(mags_b_err,chisq_r,color='purple')
plt.xlabel('Magnitude Error (r-band)')
plt.ylabel('Reduced Chi Squared (b-band)')
plt.savefig(save_dir+f'imfit_results_{fn_stub}_compare_chisqr_mag_b_err.png')
plt.clf()
plt.scatter(n_err, chisq_r,color='purple')
plt.xlabel('Sersic Index Error')
plt.ylabel('Reduced Chi Squared (r-band)')
plt.savefig(save_dir+f'imfit_results_{fn_stub}_compare_chisqr_n_err.png')
plt.clf()
plt.scatter(ell_err, chisq_r,color='purple')
plt.xlabel('Ellipticity Error')
plt.ylabel('Reduced Chi Squared (r-band)')
plt.savefig(save_dir+f'imfit_results_{fn_stub}_compare_chisqr_ell_err.png')
plt.clf()
plt.scatter(pa_err, chisq_r,color='purple')
plt.xlabel('Position Angle Error')
plt.ylabel('Reduced Chi Squared (r-band)')
plt.savefig(save_dir+f'imfit_results_{fn_stub}_compare_chisqr_pa_err.png')

if comp == 'comp_2':
    plt.clf()
    plt.scatter(radii_err, sky_r,color='purple')
    plt.xlabel('Fractional Radius Error')
    plt.ylabel('Sky Level (r-band)')
    plt.axvline(x=0, color='black', linestyle='dashed', linewidth=3)
    plt.axhline(y=0, color='black', linestyle='dashed', linewidth=3)
    plt.savefig(save_dir+f'imfit_results_{fn_stub}_compare_radius_err_I_0_r.png')
    plt.clf()
    plt.scatter(radii_err, sky_b,color='purple')
    plt.xlabel('Fractional Radius Error')
    plt.ylabel('Sky Level (b-band)')
    plt.axvline(x=0, color='black', linestyle='dashed', linewidth=3)
    plt.axhline(y=0, color='black', linestyle='dashed', linewidth=3)
    plt.savefig(save_dir+f'imfit_results_{fn_stub}_compare_radius_err_I_0_b.png')

plt.clf()
#bins = np.linspace((config['sersic_params']['r_e_min']-r_true)/r_true,(config['sersic_params']['r_e_max']-r_true)/r_true,num_bins)
plt.hist(pos_err,bins='auto',color='purple',alpha=0.8)
plt.axvline(x=0, color='black', linestyle='dashed', linewidth=3)
plt.axvline(x=np.mean(pos_err), color='firebrick', linestyle='dashed', linewidth=3, label=f'average: {round(np.mean(pos_err),3)}')
plt.xlabel('Position Error (pixels)')
plt.axvline(x=np.median(pos_err), color='gold', linestyle='dashed', linewidth=3, label=f'median: {round(np.median(pos_err),3)}')
plt.legend(prop={'size': 15})
plt.savefig(save_dir+f'imfit_results_{fn_stub}_num-iters{num_iters}_position.png')

plt.clf()
#bins = np.linspace((config['sersic_params']['r_e_min']-r_true)/r_true,(config['sersic_params']['r_e_max']-r_true)/r_true,num_bins)
plt.hist((radii-r_true)/r_true,bins='auto',color='purple',alpha=0.8)
plt.axvline(x=0, color='black', linestyle='dashed', linewidth=3)
plt.axvline(x=np.mean((radii-r_true)/r_true), color='firebrick', linestyle='dashed', linewidth=3, label=f'average: {round(np.mean((radii-r_true)/r_true),3)}')
plt.xlabel('Fractional Radius Error')
plt.axvline(x=np.median((radii-r_true)/r_true), color='gold', linestyle='dashed', linewidth=3, label=f'median: {round(np.median((radii-r_true)/r_true),3)}')
plt.legend(prop={'size': 15})
plt.savefig(save_dir+f'imfit_results_{fn_stub}_num-iters{num_iters}_radii.png')

plt.clf()
#num_bins = 20
#bins = np.linspace(config['sersic_params']['n_min']-n_true,config['sersic_params']['n_max']-n_true,num_bins)
plt.hist(n-n_true,bins='auto',color='purple',alpha=0.8)
plt.axvline(x=0, color='black', linestyle='dashed', linewidth=3)
plt.axvline(x=np.mean(n-n_true), color='firebrick', linestyle='dashed', linewidth=3, label=f'average: {round(np.mean(n-n_true),3)}')
plt.xlabel('Sersic Index Error')
plt.axvline(x=np.median(n-n_true), color='gold', linestyle='dashed', linewidth=3, label=f'median: {round(np.median(n-n_true),3)}')
plt.legend(prop={'size': 15})
plt.savefig(save_dir+f'imfit_results_{fn_stub}_num-iters{num_iters}_n.png')

plt.clf()
#num_bins = 20
#bins = np.linspace(config['sersic_params']['n_min']-n_true,config['sersic_params']['n_max']-n_true,num_bins)
plt.hist(ell-ell_true,bins='auto',color='purple',alpha=0.8)
plt.axvline(x=0, color='black', linestyle='dashed', linewidth=3)
plt.axvline(x=np.mean(ell-ell_true), color='firebrick', linestyle='dashed', linewidth=3, label=f'average: {round(np.mean(ell-ell_true),3)}')
plt.xlabel('Ellipticity Error')
plt.axvline(x=np.median(ell-ell_true), color='gold', linestyle='dashed', linewidth=3, label=f'median: {round(np.median(ell-ell_true),3)}')
plt.legend(prop={'size': 15})
plt.savefig(save_dir+f'imfit_results_{fn_stub}_num-iters{num_iters}_ell.png')

if plot_I_e:
    plt.clf()
    #num_bins = 20
    #bins = np.linspace(config['sersic_params']['I_e_min']-I_e_r_true,config['sersic_params']['I_e_max']-I_e_r_true,num_bins)
    plt.hist(I_e_r-I_e_r_true,bins='auto',color='purple',alpha=0.8)
    plt.axvline(x=0, color='black', linestyle='dashed', linewidth=3)
    plt.axvline(x=np.mean(I_e_r-I_e_r_true), color='firebrick', linestyle='dashed', linewidth=3, label=f'average: {round(np.mean(I_e_r-I_e_r_true),3)}')
    plt.xlabel('Intensity Error (r-band)')
    plt.axvline(x=np.median(I_e_r-I_e_r_true), color='gold', linestyle='dashed', linewidth=3, label=f'median: {round(np.median(I_e_r-I_e_r_true),3)}')
    plt.legend(prop={'size': 15})
    plt.savefig(save_dir+f'imfit_results_{fn_stub}_num-iters{num_iters}_I_e_r.png')

    plt.clf()
    #num_bins = 20
    #bins = np.linspace(config['sersic_params']['I_e_min']-I_e_b_true,config['sersic_params']['I_e_max']-I_e_b_true,num_bins)
    plt.hist(I_e_b-I_e_b_true,bins='auto',color='purple',alpha=0.8)
    plt.axvline(x=0, color='black', linestyle='dashed', linewidth=3)
    plt.axvline(x=np.mean(I_e_b-I_e_b_true), color='firebrick', linestyle='dashed', linewidth=3, label=f'average: {round(np.mean(I_e_b-I_e_b_true),3)}')
    plt.xlabel('Intensity Error (b-band)')
    plt.axvline(x=np.median(I_e_b-I_e_b_true), color='gold', linestyle='dashed', linewidth=3, label=f'median: {round(np.median(I_e_b-I_e_b_true),3)}')
    plt.legend(prop={'size': 15})
    plt.savefig(save_dir+f'imfit_results_{fn_stub}_num-iters{num_iters}_I_e_b.png')

plt.clf()
#num_bins = 20
#bins = np.linspace(14-mag_r_true,17-mag_r_true,num_bins)
plt.hist(mags_r-mag_r_true,bins='auto', color='purple',alpha=0.8)
plt.axvline(x=0, color='black', linestyle='dashed', linewidth=3)
plt.axvline(x=np.mean(mags_r-mag_r_true), color='firebrick', linestyle='dashed', linewidth=3, label=f'average: {round(np.mean(mags_r-mag_r_true),3)}')
plt.xlabel('Magnitude Error (r-band)')
plt.axvline(x=np.median(mags_r-mag_r_true), color='gold', linestyle='dashed', linewidth=3, label=f'median: {round(np.median(mags_r-mag_r_true),3)}')
plt.legend(prop={'size': 15})
plt.savefig(save_dir+f'imfit_results_{fn_stub}_num-iters{num_iters}_mags_r.png')

plt.clf()
#num_bins = 20
#bins = np.linspace(14-mag_b_true,17-mag_b_true,num_bins)
plt.hist(mags_b-mag_b_true,bins='auto',color='purple',alpha=0.8)
plt.axvline(x=0, color='black', linestyle='dashed', linewidth=3)
plt.axvline(x=np.mean(mags_b-mag_b_true), color='firebrick', linestyle='dashed', linewidth=3, label=f'average: {round(np.mean(mags_b-mag_b_true),3)}')
plt.xlabel('Magnitude Error (b-band)')
plt.axvline(x=np.median(mags_b-mag_b_true), color='gold', linestyle='dashed', linewidth=3, label=f'median: {round(np.median(mags_b-mag_b_true),3)}')
plt.legend(prop={'size': 15})
plt.savefig(save_dir+f'imfit_results_{fn_stub}_num-iters{num_iters}_mags_b.png')

plt.clf()
#num_bins = 20
#bins = np.linspace(0.5-color_true,2.0-color_true,num_bins)
plt.hist(colors-color_true,bins='auto',color='purple',alpha=0.8)
plt.axvline(x=0, color='black', linestyle='dashed', linewidth=3)
plt.axvline(x=np.mean(colors-color_true), color='firebrick', linestyle='dashed', linewidth=3, label=f'average: {round(np.mean(colors-color_true),3)}')
plt.xlabel('Color Error')
plt.axvline(x=np.median(colors-color_true), color='gold', linestyle='dashed', linewidth=3, label=f'median: {round(np.median(colors-color_true),3)}')
plt.legend(prop={'size': 15})
plt.savefig(save_dir+f'imfit_results_{fn_stub}_num-iters{num_iters}_colors.png')


###### SBF Magnitude stuff
if not random_sersics :
    file = open(in_dir+f'sbf_results_{fn_stub}_num-iters{num_iters}_weighted_avg_sbf_mags', "rb")
    sbf_mags = np.load(file,allow_pickle=True)
    file.close
    file = open(in_dir+f'sbf_results_{fn_stub}_num-iters{num_iters}_weighted_avg_sbf_dist_a', "rb")
    sbf_dist_a = np.load(file,allow_pickle=True)
    file.close
    file = open(in_dir+f'sbf_results_{fn_stub}_num-iters{num_iters}_weighted_avg_sbf_dist_b', "rb")
    sbf_dist_b = np.load(file,allow_pickle=True)
    file.close

    plt.clf()
    plt.scatter(sbf_mags-sbf_mag_true, chisq_r,color='purple')
    plt.xlabel('SBF magnitude Error')
    plt.ylabel('Reduced Chi Squared (r-band)')
    plt.savefig(save_dir+f'imfit_results_{fn_stub}_compare_chisqr_sbfmag_err.png')
    plt.clf()
    plt.scatter(sbf_mags-sbf_mag_true, chisq_b,color='purple')
    plt.xlabel('SBF magnitude Error')
    plt.ylabel('Reduced Chi Squared (b-band)')
    plt.savefig(save_dir+f'imfit_results_{fn_stub}_compare_chisqb_sbfmag_err.png')
    plt.clf()
    plt.scatter((sbf_dist_a/1e6)-true_dist, chisq_r,color='purple')
    plt.xlabel('Distance Error (Mpc)')
    plt.ylabel('Reduced Chi Squared (r-band)')
    plt.savefig(save_dir+f'imfit_results_{fn_stub}_compare_chisqr_dist_err.png')
    plt.clf()
    plt.scatter((sbf_dist_a/1e6)-true_dist, chisq_b,color='purple')
    plt.xlabel('Distance Error (Mpc)')
    plt.ylabel('Reduced Chi Squared (b-band)')
    plt.savefig(save_dir+f'imfit_results_{fn_stub}_compare_chisqb_dist_err.png')

    plt.clf()
    plt.hist(sbf_mags-sbf_mag_true,bins='auto',color='purple',alpha=0.8)
    plt.axvline(x=0, color='black', linestyle='dashed', linewidth=3)
    plt.axvline(x=np.mean(sbf_mags-sbf_mag_true), color='firebrick', linestyle='dashed', linewidth=3, label=f'average: {round(np.mean(sbf_mags-sbf_mag_true),3)}')
    plt.axvline(x=np.median(sbf_mags-sbf_mag_true), color='gold', linestyle='dashed', linewidth=3, label=f'median: {round(np.median(sbf_mags-sbf_mag_true),3)}')
    plt.xlabel('SBF Magnitude Error')
    plt.legend(prop={'size': 15})
    plt.savefig(save_dir+f'sbf_results_{fn_stub}_num-iters{num_iters}_sbfmagnitude.png')

    plt.clf()
    plt.hist((sbf_dist_a/1e6)-true_dist,bins='auto',color='purple',alpha=0.8,range=(-1,2))
    plt.axvline(x=0, color='black', linestyle='dashed', linewidth=3)
    plt.axvline(x=np.mean((sbf_dist_a/1e6)-true_dist), color='firebrick', linestyle='dashed', linewidth=3, label=f'average: {round(np.mean((sbf_dist_a/1e6)-true_dist),3)}')
    plt.axvline(x=np.median((sbf_dist_a/1e6)-true_dist), color='gold', linestyle='dashed', linewidth=3, label=f'median: {round(np.median((sbf_dist_a/1e6)-true_dist),3)}')
    plt.xlabel('SBF Distance Error A (Jerjen Eqn. 2)')
    plt.legend(prop={'size': 15})
    plt.savefig(save_dir+f'sbf_results_{fn_stub}_num-iters{num_iters}_distance_a.png')

    plt.clf()
    plt.hist((sbf_dist_b/1e6)-true_dist,bins='auto',color='purple',alpha=0.8, range=(-1,2))
    plt.axvline(x=0, color='black', linestyle='dashed', linewidth=3)
    plt.axvline(x=np.mean((sbf_dist_b/1e6)-true_dist), color='firebrick', linestyle='dashed', linewidth=3, label=f'average: {round(np.mean((sbf_dist_b/1e6)-true_dist),3)}')
    plt.axvline(x=np.median((sbf_dist_b/1e6)-true_dist), color='gold', linestyle='dashed', linewidth=3, label=f'median: {round(np.median((sbf_dist_b/1e6)-true_dist),3)}')
    plt.xlabel('SBF Distance Error B (Jerjen Eqn. 3)')
    plt.legend(prop={'size': 15})
    plt.savefig(save_dir+f'sbf_results_{fn_stub}_num-iters{num_iters}_distance_b.png')


plt.clf()
plt.scatter(mags_r-mag_r_true, sbf_mags-sbf_mag_true,color='purple')
plt.xlabel('Magnitude Error (r-band)')
plt.ylabel('SBF Magnitude Error')
plt.savefig(save_dir+f'sbf_results_{fn_stub}_mag_r_err_sbfmag_err.png')

plt.clf()
plt.scatter(mags_b-mag_b_true, sbf_mags-sbf_mag_true,color='purple')
plt.xlabel('Magnitude Error (b-band)')
plt.ylabel('SBF Magnitude Error')
plt.savefig(save_dir+f'sbf_results_{fn_stub}_mag_b_err_sbfmag_err.png')

plt.clf()
plt.scatter(colors-color_true, sbf_mags-sbf_mag_true,color='purple')
plt.xlabel('Color Error')
plt.ylabel('SBF Magnitude Error')
plt.savefig(save_dir+f'sbf_results_{fn_stub}_color_err_sbfmag_err.png')

plt.clf()
plt.scatter(colors-color_true, (sbf_dist_a/1e6)-true_dist,color='purple')
plt.xlabel('Color Error')
plt.ylabel('SBF Distance Error A (Jerjen Eqn. 2)')
plt.savefig(save_dir+f'sbf_results_{fn_stub}_color_err_dist_err.png')

plt.clf()
plt.scatter(mags_r-mag_r_true, (sbf_dist_a/1e6)-true_dist,color='purple')
plt.xlabel('Magnitude Error (r-band)')
plt.ylabel('SBF Distance Error A (Jerjen Eqn. 2)')
plt.savefig(save_dir+f'sbf_results_{fn_stub}_mag_r_err_dist_err.png')


'''
bad_rad_mask = abs(radii_err) > 10.
good_rad_mask = ~bad_rad_mask
really_good_rad_mask = abs(radii_err) < 2.
really_bad_large_rad_mask = radii_err > 30.
really_bad_small_rad_mask = radii_err < -30.

print('BAD RADIUS MEASUREMENT: ',iters[bad_rad_mask])
print('GOOD RADIUS MEASUREMENT: ',iters[good_rad_mask])
print('REALLY GOOD RADIUS MEASUREMENT: ',iters[really_good_rad_mask])
print('REALLY BAD (Large) RADIUS MEASUREMENT: ',iters[really_bad_large_rad_mask])
print('REALLY BAD (Small) RADIUS MEASUREMENT: ',iters[really_bad_small_rad_mask])
'''
