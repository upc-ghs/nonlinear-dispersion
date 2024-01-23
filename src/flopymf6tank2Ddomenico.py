'''
FloPy configuration for the model of transport through 
a two-dimensional homogeneous porous medium. Dispersion 
model employs the expression in Chiogna et al. (2010) 
which is given to MODFLOW by providing an equivalent 
dispersivity.

@references
    Chiogna et al. (2010), Evidence of Compound-Dependent Hydrodynamic and 
    Mechanical Transverse Dispersion by Multitracer Laboratory Experiments
'''

# Dependencies
import os
import sys
import numpy as np
import flopy
import matplotlib.pyplot as plt
import argparse
import xmipy
import scipy.special as spc

# self
import config
import dispersion

#############
# Arguments #
#############
parser = argparse.ArgumentParser( description='Setup simulations for experiment.' )
parser.add_argument( '--write'    , action='store_true', help='write simulation files.' )
parser.add_argument( '--run'      , action='store_true', help='executes simulation, requires mf6lib.' )
parser.add_argument( '--force'    , action='store_true', help='skip the prompt after verification of time step.' )
parser.add_argument( '--cprofile' , action='store_true', help='displays a transverse profile of concentration.' )
args   = parser.parse_args()


#########################
# General configuration #
#########################
HIGHRES        = False # adopt high resolution parameters
targetvelocity = 10 # m/day 

# total time in hours 
if targetvelocity == 0.5: 
    tlength       = 20.0   
if targetvelocity == 1: 
    tlength       = 10.0   
if targetvelocity == 10: 
    tlength       = 1.0   
if targetvelocity == 20: 
    tlength       = 0.5

if HIGHRES:
    tn_time_steps = 15000
else:
    tn_time_steps = 1500


stressperiods = [
        {
            'id': 0,
            'length'       : tlength,
            'n_time_steps' : tn_time_steps, 
            'ts_multiplier': 1,
            'steady_state' : True,
            'printfreqts'  : 1, 
        },
    ]
experiment_folder = config.FOLDERS['sims']
base_model_name   = 'mf6model'

# component
injconc           = 1
dpdsdiff          = 0.1 # dp/ds
sim_name          = 'NLDOM2D'
davg              = 1e-9*1e4*3600   # avg diffusion, the same for all sims
porosity          = 0.35
delta             = 6
betaT             = 0.5
ck                = 0.01             # np.sqrt(ms) ( Hazen, 1892 )
diam              = 0.0015           # meters (1.5mm)
hk                = diam**2/(ck**2)  # ( Hazen, 1892 )
hk                = hk*100*3600      # m/s
sim_ws            = os.path.join( experiment_folder, sim_name )

print( ' - INFO:  Grain diameter : ', diam*1000  , 'mm   ' ) 
print( ' - INFO:  Hydraulic cond : ', hk         , 'm/s  ' )
print( ' - INFO:  Hydraulic cond : ', hk         , 'cm/hr' )


# Init sim and configure stress periods
sim = flopy.mf6.MFSimulation(
        sim_name=sim_name,
        sim_ws=sim_ws,
        exe_name='mf6'
    )
# tdis package
perioddata = []
for sp in stressperiods:
    perioddata.append( [sp['length'], sp['n_time_steps'], sp['ts_multiplier']] )
tdis = flopy.mf6.ModflowTdis(
    sim,
    nper=len(stressperiods),
    perioddata=perioddata,
    time_units='hours',
)

# IMS
hclose = 1e-5
ninner = 100
nouter = 50
relax  = 0.97
complexity = 'COMPLEX'
ims = flopy.mf6.ModflowIms(
        sim,
        pname               = 'ims',
        print_option        = 'SUMMARY',
        complexity          = complexity,
        outer_dvclose       = hclose, 
        outer_maximum       = nouter, 
        inner_maximum       = ninner, 
        inner_dvclose       = hclose, 
        linear_acceleration = 'BICGSTAB',
    )

# gwf
gwf = flopy.mf6.ModflowGwf(
        sim,
        modelname=base_model_name,
        save_flows=True
    )
# dis package 
LTANK  = 40.0  # cm
nlay   = 1
if HIGHRES:
    wcells = 31    # number of cells injecting
    nrow   = 250
    ncol   = 1250
else:
    wcells = 12    # number of cells injecting
    nrow   = 160
    ncol   = 800
delc   = LTANK/ncol
delz   = delc
delr   = delc
dis    = flopy.mf6.ModflowGwfdis(
            gwf,
            nlay=nlay, nrow=nrow, ncol=ncol,
            delr=delr, delc=delc, top=delz, botm=0,
            length_units='CENTIMETERS',
        )

# npf package
npf = flopy.mf6.ModflowGwfnpf(
    gwf,
    save_specific_discharge=True,
    icelltype=0,
    k=hk,
)

# sto package
sto = flopy.mf6.ModflowGwfsto(gwf, ss=0,sy=0, steady_state=True)

# chd 
# Darcy's law
print( ' - INFO:  Target velocity : ', targetvelocity, '  m/d  ')
targetvelocity  = targetvelocity*100/24   # cm/hr
print( ' - INFO:  Target velocity : ', targetvelocity, ' cm/hr ')
qflow           = targetvelocity*porosity # cm/hr
gradh           = qflow/hk                # deltah/length
deltah          = gradh*LTANK
href            = 10 # could be anything 
hdin            = href + deltah
hdout           = href
totalflowrate   = qflow*delc**2*nrow
print( ' - INFO:  Head gradient  : ', gradh  , ' -  ')
print( ' - INFO:  Delta head     : ', deltah , ' cm ')
print( ' - INFO:  Head inlet     : ', hdin   , ' cm ')
print( ' - INFO:  Head outlet    : ', hdout  , ' cm ')
print( ' - INFO:  Flow rate      : ', totalflowrate , ' cm3/hr ' )
print( ' - INFO:  Pore volume    : ', porosity*nrow*ncol*nlay*delc**3/totalflowrate    , ' hr  ')
print( ' - INFO:  Pore volume    : ', porosity*nrow*ncol*nlay*delc**3/totalflowrate*60 , ' min ')

# assign chd for in flow cells 
# center cells inject secondary species
# the rest inject the primary species
chdspd     = []
thesecells = [ int(0.5*(nrow - wcells)) + iw for iw in range(wcells) ]
for ir in range(nrow):
    if ir in thesecells:
        chdspd.append( [ (0, ir, 0), hdin, injconc ] ) # inject secondary at the center
    else:
        chdspd.append( [ (0, ir, 0), hdin, 0       ] )     # inject primary in the rest
chdin = flopy.mf6.ModflowGwfchd(
        gwf,
        stress_period_data=chdspd,
        maxbound=len(chdspd),
        auxiliary="CONCENTRATION",
        pname='CHD-1'
    )
# assign chd for outflow cells 
chdspd     = []
for ir in range(nrow):
    chdspd.append( [ (0, ir, ncol-1), hdout ] )
chdout = flopy.mf6.ModflowGwfchd(
        gwf,
        stress_period_data=chdspd,
        maxbound=len(chdspd),
        pname='CHD-2'
    )

# ic (linear gradient)
icst = np.zeros(shape=(nlay,nrow,ncol))
xarr = np.arange(0.5*delc,LTANK+0.5*delc,delc)
for ir in range(nrow):
    icst[0,ir,:] = -1*gradh*xarr + hdin
flopy.mf6.ModflowGwfic(
        gwf,
        strt=icst,
    )

# oc
head_file   = base_model_name + '.hds'
budget_file = base_model_name + '.bud'   
nsteps      = stressperiods[0]['n_time_steps']
tlength     = stressperiods[0]['length']
deltat      = tlength/nsteps
printfreqts = stressperiods[0]['printfreqts']
times       = np.arange(0,tlength+deltat,deltat)
printtimes  = times[::printfreqts]
saverecord  = [ ('HEAD','LAST'), ('BUDGET','LAST')]
flopy.mf6.ModflowGwfoc(
    gwf,
    head_filerecord  =head_file,
    budget_filerecord=budget_file,
    saverecord=saverecord
)

#######
# GWT #
#######
gwtname = "gwt-" + base_model_name
gwt     = flopy.mf6.MFModel(
        sim,
        model_type="gwt6",
        modelname=gwtname,
        model_nam_file="{}.nam".format(gwtname),
    )
gwt.name_file.save_flows = True

# create iterative model solution and register the gwt model with it
imsgwt = flopy.mf6.ModflowIms(
    sim,
    print_option        = "SUMMARY",
    outer_dvclose       = hclose,
    outer_maximum       = nouter,
    under_relaxation    = "NONE",
    inner_maximum       = ninner,
    inner_dvclose       = hclose,
    linear_acceleration = "BICGSTAB",
    complexity          = complexity,
    filename            = "{}.ims".format(gwtname),
)
sim.register_ims_package(imsgwt, [gwt.name])

# disgwt
disgwt = flopy.mf6.ModflowGwfdis(
        gwt,
        nlay=nlay, nrow=nrow, ncol=ncol,
        delr=delr, delc=delc, top=delz, botm=0,
        length_units='CENTIMETERS',
    )

# adv
scheme = "TVD"
flopy.mf6.ModflowGwtadv(
        gwt,
        scheme=scheme,
        filename="{}.adv".format(gwtname)
    )

#ic conc
#saturated with primary species
sconc = 0
flopy.mf6.ModflowGwtic(
    gwt,
    strt=sconc,
    filename="{}.ic".format(gwtname)
)

# some parameters for dispersion
deltax      = delc 
tscalediff  = deltax**2/(davg)    # cm2/hr
tscaleadv   = deltax/targetvelocity  # hr 
ct          = 0.1
displ       = 0.5*diam*100*targetvelocity                 # cm2/hr
displphi    = porosity*davg + 0.5*diam*100*targetvelocity # cm2/hr
tscaledl    = deltax**2/(displphi)                        # cm2/hr
trec        = np.min([ct*np.mean(tscaledl), ct*tscaleadv])
peclet      = targetvelocity*diam*100/davg

# dispersion parameters
alphal  = 0.5*diam*100 # cm
alphat  = (diam*100/peclet)*dispersion.nonlinear_dispersion_function( peclet, delta, betaT )
dmeff   = porosity*davg

# dsp package
flopy.mf6.ModflowGwtdsp(
    gwt,
    alh      = alphal,
    ath1     = alphat,
    diffc    = dmeff, 
    xt3d_off = False,
    filename = "{}.dsp".format(gwtname),
)


print( ' # Dispersion analysis:    ' )
print( '   - deltat (given)      : ', stressperiods[0]['length']/stressperiods[0]['n_time_steps'] )
print( '   - deltat (recommended): ', trec )
print( '   - nsteps (given)      : ', stressperiods[0]['n_time_steps'] ) 
print( '   - nsteps (recommended): ', stressperiods[0]['length']/trec ) 
print( '   - effective CT        : ', (stressperiods[0]['length']/stressperiods[0]['n_time_steps'])/trec*ct  ) 
print( 'Continue ... ? (press c to continue, or q to exit)' )
if not args.force:
    import pdb
    pdb.set_trace()


# mst package
flopy.mf6.ModflowGwtmst(
    gwt,
    porosity          = porosity,
    first_order_decay = False,
    decay             = None,
    decay_sorbed      = None,
    sorption          = None,
    bulk_density      = None,
    distcoef          = None,
    filename="{}.mst".format(gwtname),
)

# ssm
sourcerecarray = [('CHD-1','AUX','CONCENTRATION')]
flopy.mf6.ModflowGwtssm(
    gwt,
    sources=sourcerecarray,
    filename="{}.ssm".format(gwtname),
)

# gwt oc 
concentration_file = "{}.ucn".format(gwtname)
saverecord         = [ ('CONCENTRATION', 'FIRST'), ('CONCENTRATION','LAST') ]
flopy.mf6.ModflowGwtoc(
    gwt,
    concentration_filerecord=concentration_file,
    concentrationprintrecord=[
        ("COLUMNS", 10, "WIDTH", 15, "DIGITS", 6, "GENERAL")
    ],
    saverecord=saverecord, 
    filename="{}.oc".format(gwtname),
)

# gwfgwt
flopy.mf6.ModflowGwfgwt(
    sim,
    exgtype="GWF6-GWT6",
    exgmnamea=base_model_name,
    exgmnameb=gwtname,
    filename="{}.gwfgwt".format(base_model_name),
)

#############
# Write MF6 #
#############
# No need to use the MODFLOW API for this problem
#################################################
write_simulation = args.write
if write_simulation:
    sim.write_simulation()
if args.run: 
    sim.run_simulation()


#############
# C profile #
#############
if args.cprofile:
    '''
    Generates a figure of the profile of concentration 
    compared to analytical solution of Domenico. 

    The figure is stored as png in the folder of the simulation (sim_ws)
    '''

    # Dependencies
    import matplotlib.colors as colors
    import matplotlib.cm as cm 
    # self
    import analytical

    gwtname  = gwtname.lower()
    layer_id = 0
    fname    = os.path.join( os.getcwd(), sim_ws, gwtname+'.ucn')
    ucnobj   = flopy.utils.HeadFile(fname, text='CONCENTRATION', precision='double')
    data     = ucnobj.get_alldata()
    xindex   = -100 

    # Figure handling
    fig            = plt.figure()
    modelgrid      = gwf.modelgrid
    xposition      = modelgrid.xycenters[0][xindex]
    domainwidth    = nrow*delr
    npoints        = nrow
    deltay         = delr
    yvector        = modelgrid.xycenters[1] - domainwidth*0.5
    injectionwidth = wcells*deltay 
    velocity       = targetvelocity
    dispersiony    = dmeff + alphat*velocity
    solution       = analytical.domenico(
            xposition,
            yvector,
            injectionwidth,
            dispersiony,
            velocity,
            bounded=True,
            domainwidth=domainwidth,
        )
    plt.plot( yvector, solution, color='b', label='analytical' )
    plt.plot( yvector, data[-1,0,:,xindex], color='r', label='numerical' )
    plt.legend()
    ax = plt.gca()
    ax.set_xlabel('y')
    ax.set_ylabel('C')
    plt.savefig( os.path.join( sim_ws, 'cprof'+sim_name+'.png') )
