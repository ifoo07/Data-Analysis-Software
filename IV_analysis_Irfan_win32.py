# -* coding UTF-8 *-
'''
Change: deals with nlog files (recalculated)

DESCRIPT: New code for AFMD solar cell analysis. Cleaned up/improved usage.
Generates Box plots for each experimental variable (bin).
Works with several log files (2-1 onwards).

USER instructions:
- Put all log files you care about in a folder
- In this folder add a bins csv file, eg myBins.csv, 
    following the template (myBinsEg.csv)
- Follow Instructions (make sure bins.csv is ok)

NB:
- From ECHO you can have two different stripes on the same quad. 
    These correspond to odd and even pixels in the snaith meas scheme 
- ALWAYS ALIGN + (odd) to proper corner (right side looking from computer) in solar sim.
- The cull duplicates option will take the highest eff in case a pixel has been measured twice
Alternative:

Programmer notes:
1. Load settings (myBins.csv)
2. Generate panda data base
3. Make Bin column
4. Group and Box plot

Author: I Ramirez, AFMD 01/2016
Changes: graphs between 0 and 1.05*max(value)
Edit: Anna's version edited by Irfan to run on win32 with updated pandas (ver 1.3)

Runs python 3 pyqt4
'''

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import os, sys
import pandas as pd
import math
import time
import seaborn as sns

import warnings
warnings.filterwarnings("ignore")

#rescale_area = 0.0919/0.01


#user parameters
autoAssignBin = False
plot_best_IVs = True
physics_user = 'jungbluth'

plot_light_of_best = True
plot_dark_of_best = False
flip_iv = True # for p-i-n vs n-i-p
logplot = False
boxplots = False

#internal parameters
logStruct = ['Jsc','Eff', 'Voc', 'FF', '?', 'sCurr','sEff','sRatio', 'label', 'var', 'value', 'pix', 'Intensity', '??', 'path']
columns = ['quad', 'pix', 'logF', 'exp_var', 'Voc', 'Jsc', "FF", 'Eff', 'scan', 'bin' ,'path'] # things we care about

df = pd.DataFrame(index= [], columns=columns) #empty data frame :)
saved_state = 0, #because we have to worry about scan direction,

# CHANGE THIS!
#mount = '/Volumes/CMGroups/HJSGroup'
mount = '/home/jungbluth/Ydrive/Groups/Condensed Matter/HJSGroup'
#mount_local = '~/Desktop/ZnPc_Data/ZnPc_series_JV_data'
mount_local = 'C:/Users/Irfan Habib/Desktop/JV data'

def parse(f):
    '''load .log .nlog and bins properly'''
    if f.endswith('.nlog'): #get rid of extra columns that would mess stuff up 
        return pd.read_csv(f, sep = '\t', index_col= None) #sep = '\t'

    elif f.endswith('.log'):
        try:
            dtemp =  pd.read_csv(f, sep = '\t', header = None, index_col= None)
            dtemp.columns = logStruct #important, add index here
        except:
            dtemp =  pd.read_csv(f, sep = ',', header = None, index_col= None)
            dtemp.columns = logStruct #important, add index here
        return dtemp
    else:
        return pd.read_csv(f)

def processRow(row, logF):
    '''add row to main data frame'''
    global df, setsDf
    sDf = setsDf[setsDf.logFileAlias == logF] #the frame appropriate for this logF
    #this is required so code works accross different log files with potential label repeats
    #add stuff in order to list which we append to df. --> list order !

    if row['path'].endswith('div1'): return #if dark scan skip
    if row['path'].endswith('div2'): return #if dark scan skip
    # print ('row:',row)
    quad = row['label']
    pix = row['pix']
    l = [quad, pix]
    l.append(logF)
    l.append(row['var'])
    scanValues = row['Voc'],row['Jsc'], row['FF'], row['Eff']
    l+= [abs(float(r)) for r in scanValues]
    l.append(assertScan(quad, pix, row)) #scan dir

    # if row['path'].endswith('div1') or row['path'].endswith('div2'):
    #     return #if dark scan skip
    # else:
    #     quad = row['label']
    #     pix = row['pix']
    #     l = [quad, pix]
    #     l.append(logF)
    #     l.append(row['var'])
    #     scanValues = row['Voc'],row['Jsc'], row['FF'], row['Eff']
    #     l+= [abs(float(r)) for r in scanValues]
    #     l.append(assertScan(quad, pix, row)) #scan dir

    #work out what the correct bin is from loaded bins file
    if autoAssignBin:
        l.append(quad)

    else:
        parity = int(pix)%2 #which stripe
        if parity == 0: parity = 'even'
        else: parity = 'odd'
        sDf_ = sDf[sDf.Stripe == parity].copy()
        bin = sDf_[sDf_.quad == quad]['variable'].to_string(index=False)#...?
        print('bin:', bin)
        l.append(bin)
        l.append(row['path'])
        # try: #assign bin by comparing to desired stuff
        #     bin = sDf[sDf.Stripe == parity].loc[quad]['variable']#...?
            
        #     l.append(bin)
        #     l.append(row['path'])
        # except Exception as e: print(e) 
        # # KeyError: pass #not found
    # print('l:',l,len(l))
    df_temp = pd.DataFrame(dict(zip(columns, l)), index =range(1)) #make into frame
    print(l)
    df =df.append(df_temp, ignore_index=True) #add to main frame
def assertScan(q, p, row):
    '''function to work out scan type'''
    global saved_state
    qp = str(q)+str(p)
    if float(row['sEff']) != 0 and float(row['sCurr']) != 0:
        saved_state = qp
        return 'stab'
    elif qp == saved_state: #if previous was already same pixel then this is back scan 
        saved_state = qp
        return 'second'
    else:
        saved_state = qp
        return 'first'

def boxplot(datF, field):
    '''http://stackoverflow.com/questions/23519135/dot-boxplots-from-dataframes'''
    datF.boxplot(grid = False, by = 'bin', column = field)
    plt.tight_layout()
    plt.suptitle("")#remove the auto title
    bins = pd.unique(datF.sort_values(by='bin')['bin'].ravel())
    colors = plt.cm.winter(np.linspace(0,1,len(bins)))
    for i, b in enumerate(bins):
        y = np.array(datF[datF.bin == b][field])
        x = np.random.normal(i+1, 0.04, len(y))
        if field == 'FF':
            maxi = 1 # sets the y axis limit to 1 for FF

        # CHANGE THIS IF YOU WANT TO LIMIT Y-AXIS
        elif field == 'Jsc':
            maxi = 5

        else:
            maxi= max(1.05*df[field]) # sets the maximal value for Eff, Jsc and Voc (set back to 1.05?)
        plt.ylim(0.0,maxi) # set y axis limits

        plt.plot(x, y, mfc = colors[i], mec='k', ms=7, marker="o", linestyle="None")

def pruneData(datF, filename):
    """
    Function to select and discard data
    :param datF: dataFrame of JV scans [dataFrame]
    :param filename: name of (bins/pruned) file [str]
    :return:
    """

    use_file = input('Start from file? (0/1)\n>')

    bins_ = pd.unique(datF['bin'])
    # print('datF:',datF)
    # print('bins_:', bins_)
    # print(datF['bin'].ravel())
    bins = []
    for b in bins_: # This is to remove Nans from bins_
        if type(b) == str:
            bins.append(b)
    print(bins)
    if int(use_file): # Ask for file name / directory if starting from file
        binFile, dir_ = getUserSettings(no=1)
        os.chdir(dir_)
        pruned_df = pd.read_csv(binFile, sep=',', header=0)
        # print(pruned_df)

    else: # Plot and sort JV scans if not starting from file

        if filename.endswith('.csv'):
            filename = filename.split('.')[0]

        # bins = pd.unique(datF.sort_values(by='bin')['bin'].ravel())
        print('bins:', bins)
        pruned_df = pd.DataFrame()
        # print(datF)
        for i, b in enumerate(bins): # iterate through each sample
            y_df = datF[datF.bin == b].copy()
            y_df = y_df[df['pix']!=3].copy() # remove small pixels
            y_df = y_df[df['pix']!=6].copy()
            
            keep_list = []

            for iv_file in y_df['path']:
                iv_path = fix_path(iv_file)
                
                if iv_path:
                    # Load data
                    v, j = loadxIV(iv_path)

                    iv_path.replace("\\","/")

                    label_ = iv_path.split("/")[-1]
                    ls = '-'

                    # Plot data
                    plt.figure()
                    plt.tight_layout()
                    plt.axhline(0, color='xkcd:black', linewidth=2)
                    plt.minorticks_on()
                    plt.title(b, fontsize=15)
                    plt.xlabel('Voltage / V', fontsize=15, fontweight='medium')
                    plt.ylabel('Current Density / mA/$\mathregular{cm^2}$', fontsize=15, fontweight='medium')

                    plt.grid(False)
                    plt.rcParams['figure.facecolor'] = 'xkcd:white'
                    plt.rcParams['figure.edgecolor'] = 'xkcd:white'
                    plt.tick_params(labelsize=15, direction='in', axis='both', which='major', length=8, width=2)
                    plt.tick_params(labelsize=15, direction='in', axis='both', which='minor', length=4, width=2)

                    plt.plot(v, j, label=label_, linewidth=5, linestyle=ls)

                    # plt.ylim(-10, 10)
                    # plt.xlim(-1, 1.7)
                    plt.legend(loc='best', prop={'size': 10})
                    plt.show(block=False)

                    # CHANGE THIS TO CHANGE TIMEFRAME TO CLOSE WINDOW
                    plt.pause(2)
                    plt.close()

                    keep = input('Keep scan? (1/0): ')
                    if keep.isnumeric():
                        keep_list.append(int(keep))

            # keep_list = [0,1,0,0,1,1,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1] # Used for testing purposes

            y_df['Use'] = keep_list

            pruned_df = pruned_df.append(y_df)

        # keep_list = [0, 0, 1, 0, 0, 1, 1, 0, 1, 1, 1, 1, 1, 0, 0, 0, 1, 1, 0, 1, 1, 1, 1, 1]
        # pruned_df['Use'] = keep_list # Used for testing purposes

        pruned_df.to_csv('~/Desktop/%s_stats.csv' % filename)
    
    select_df = pruned_df[pruned_df['Use']==1].copy() # Compile new dataframe with selected data

    bins = pd.unique(select_df['bin'].ravel())

    for u, param in enumerate(['Jsc', 'Voc', 'Eff', 'FF']):

        print("  ")
        print('------------------------------------------------')
        print('%s ' % param,str(units_terminal[u]))
        print('------------------------------------------------')

        fig, ax = plt.subplots()

        # select_df.boxplot(grid=False, by='bin', column=param)
        plt.tight_layout()
        plt.title("")
        plt.suptitle("")
        plt.rcParams['figure.facecolor'] = 'xkcd:white'
        plt.rcParams['figure.edgecolor'] = 'xkcd:white'
        ax.tick_params(labelsize=15, direction='in', axis='both', which='major', length=8, width=2)
        ax.tick_params(labelsize=15, direction='in', axis='both', which='minor', length=4, width=2)
        #colors = plt.cm.plasma(np.linspace(0, 1, len(bins)))

        for i, b in enumerate(bins):

            df_param = select_df[select_df.bin == b].copy()
            y = np.array(df_param[param])
            x = np.random.normal(i, 0.04, len(y))

            if param == 'FF':
                plt.ylabel('Fill Factor', fontsize=15, fontweight='medium')
                ax.set_ylim(0,1)

            elif param == 'Voc':
                plt.ylabel('Open-Circuit Voltage / V', fontsize=15, fontweight='medium')
            elif param == 'Jsc':
                plt.ylabel('Short-Circuit Current Density / mA/$\mathregular{cm^2}$', fontsize=15, fontweight='medium')
            elif param == 'Eff':
                plt.ylabel('Efficiency', fontsize=15, fontweight='medium')

            plt.xlabel("")
            ax.boxplot(y, positions=[i])
            ax.plot(x, y, mfc='navy', mec='k', ms=9, marker="o", linestyle="None")

            print('Statistics for bin %s : ' % b)
            print('------------------------------')
            print('Median : ', np.round(np.percentile(y, 50, interpolation='midpoint'),3))

            mean, stdv, outlier = remove_outliers(y)
            print('Mean : ', np.round(np.mean(y),3), ' +/- ', np.round(np.std(y),3))
            print('Mean (Outliers removed) : ', np.round(mean,3), ' +/- ', np.round(stdv,3))
            print('------------------------------')
            print("  ")

        ax.set_xticks(range(i + 1))
        ax.set_xticklabels(bins, fontsize=12)

    plt.figure()
    plt.tight_layout()
    plt.axhline(0, color='xkcd:black', linewidth=2)
    plt.minorticks_on()
    plt.xlabel('Voltage / V', fontsize=15, fontweight='medium')
    plt.ylabel('Current Density / mA/$\mathregular{cm^2}$', fontsize=15, fontweight='medium')

    plt.grid(False)
    plt.rcParams['figure.facecolor'] = 'xkcd:white'
    plt.rcParams['figure.edgecolor'] = 'xkcd:white'
    plt.tick_params(labelsize=15, direction='in', axis='both', which='major', length=8, width=2)
    plt.tick_params(labelsize=15, direction='in', axis='both', which='minor', length=4, width=2)

    for x, b in enumerate(bins):

        df_param = select_df[select_df.bin == b].copy()

        print('No. scans for sample %s : ' % b, str(len(df_param)))

        samples = [x[:-5] for x in df_param['path']]
        unique_samples = set(samples)

        print('No. unique devices for sample %s : ' % b, len(unique_samples))

        # CHANGE THIS TO ADD MORE COLORS FOR SCANS
        colors = ['#000000', '#801355', '#E6740C', '#F1C211', '#52A5C5']

        for i, iv_file in enumerate(df_param['path']):
            iv_path = fix_path(iv_file)

            if iv_path:
                v, j = loadxIV(iv_path)
                plt.plot(v, j, color=colors[x], linewidth=6, alpha=0.2)

    print('------------------------------')
    inds = select_df.groupby(['bin'])['Eff'].transform(max) == select_df['Eff']
    champs = select_df[inds]

    for n, iv_file in enumerate(champs['path']):
        iv_path = fix_path(iv_file)
        

        if iv_path:
            v, j = loadxIV(iv_path)
            label_ = champs[champs.path == iv_file]['bin'].values[0]
            if label_ == '$\mathregular{C_{60}}$':
                ls = '--'
            else:
                ls = '-'
            plt.plot(v, j, label=label_, linewidth=4, color=colors[n], ls=ls)

    plt.xlim(0, 1.5)
    plt.ylim(-6, 10)

    plt.legend(loc='best',prop={'size':15})
    plt.show()

def remove_outliers(data):
    """
    :param data: array or list of data
    :return: mean: mean of data without outliers
             stdv: standard deviation of data without outliers

    """

    Q1 = np.percentile(data, 25, interpolation='midpoint')
    Q3 = np.percentile(data, 75, interpolation='midpoint')
    IQR = Q3 - Q1

    min_ = Q1 - 1.5*IQR
    max_ = Q3 + 1.5*IQR

    out_removed = []
    out = []

    for value in data:
        if min_<= value <= max_:
            out_removed.append(value)
        else:
            out.append(value)

    mean = np.mean(out_removed)
    stdv = np.std(out_removed)

    return mean, stdv, out

def descriptOut():
    '''Generates a text file of the settings used and where the log files are stored'''
    pass

def cullDuplicates(df):
    '''function to remove all duplicates (keeps best eff at same time)'''
    Eff_max = df.groupby(['quad', 'pix', 'scan', 'logF']).Eff.transform(max) 
    #above: series in which duplicate effs have been replaced with maximum (still same n rows as df)
    # from http://stackoverflow.com/questions/32093829/pythonpandas-removing-duplicates-based-on-two-columns-keeping-row-with-max-va    
    return df[df.Eff == Eff_max].drop_duplicates() #drop as above only removes different measurements of same pix
    # or df.sort('C').drop_duplicates(subset=['A', 'B'], take_last=True)

def loadxIV(f):
    '''function to load IV data from liv or div files'''
    v, j = np.genfromtxt(f, skip_footer = 12).T
    if flip_iv: v*= -1

    return v,j

def fix_path(p):
    '''so doesnt just work with networked computers'''

    if sys.platform == 'linux':
        # CHANGE THIS TO MATCH YOUR PATH
        p = p.replace('\\', '/')
#        p = p.replace('H:','/Volumes/%s'%physics_user)
#        p = p.replace('/physics.ox.ac.uk/DFS', '/Volumes/winfe.physics.ox.ac.uk') #winfe
#        p = p.replace('AnnaJungbluth', 'annajungbluth') #uncomment this one when moving to sshfs ???
        p = p.replace('//groupfs.physics.ox.ac.uk/CMGroups/HJSGroup', mount)
#        p = p.replace('//groupfs.physics.ox.ac.uk/CMGroups/HJSGroup/Solar Simulator Data/annajungbluth', mount_local)
#        p = p.replace('C:', mount)
        p = p.replace('C:/Solar Simulator Data/annajungbluth', mount_local)
        p = p.replace('C:/Solar Simulator Data/AnnaJungbluth', mount_local)
    elif sys.platform =='win32':
        p = p.replace('\\', '/')
        p = p.replace('C:/Solar Simulator Data/Irfan Habib', mount_local)
        p = p.replace('C:\Solar Simulator Data\Irfan Habib', mount_local)

    if os.access(p, os.R_OK):
        if os.path.exists(p):  
            return p
        else:
            print('Error trying to plot best IVs: Path does not exist: %s'%p)
    else: 
        print('Error in trying to plot best IVs: no access to path: %s'%p)

    return None


def plot_best_cells(df_with_bins):
    '''find and plot the most efficient cell for each bin (user defined var)'''
    df_ = df_with_bins
    inds = df_.groupby(['bin'])['Eff'].transform(max) == df_['Eff']
    champs = df_[inds]

#    **** Colors ****

    # CHANGE THIS TO SET DIFFERENT COLORS
    
#    colors = ['#ffcc66', '#feb649', '#ff9e2f', '#ff8416', '#ff6600', '#A9A9A9'] #Orange
#    colors = ['#ff3300', '#df0427', '#b90034', '#900038', '#660033', '#A9A9A9'] #Red
#    colors = ['#99ccff', '#6699dd', '#3b67b9', '#1a3692', '#000066', '#A9A9A9'] #Blue
#    colors = ['#66ccff', '#10b7ff', '#009fff', '#0085ff', '#0066ff', '#A9A9A9'] #Light blue
#    colors = ['#3366cc', '#264db3', '#1a3599', '#0e1c80', '#000066', '#A9A9A9'] #Dark blue
#    colors = ['#00cc66', '#10a249', '#0e7b2f', '#055617', '#003300', '#A9A9A9'] #Green
#    colors = ['#00ff00', '#19f21a', '#25e525', '#2dd92d', '#33cc33', '#A9A9A9'] #Light green
#    colors = ['#009933', '#037e25', '#026419', '#004b0d', '#003300', '#A9A9A9'] #Dark green

#    colors = ['#feb649', '#df0427', '#ff8416', '#900038']
#    colors = ['#000000', '#660033', '#ffcc66', '#df0427', '#ff8416']
#    colors = ['#000000',  '#0710d6', '#9800ba', '#6f00c8', '#c91ba3']
#    colors = ['#000000', '#ffc700', '#f8a000', '#ee7800', '#ff3300', '#cd0808']
#    colors = ['#000000', '#45818eff', '#6aa84fff', '#f1c232ff', '#e69138ff']

    colors = ['#000000', '#F1C211', '#E6740C', '#801355', '#52A5C5']

    colors = ['#000000', \
              '#ffcc66', \
              '#ff6600', \
              '#ff9e2f', \
              '#df0427', \
              '#550091ff', \
              '#910084ff', \
              '#910036ff', \
              '#99ccff', \
              '#6699dd', \
              '#3b67b9', \
              '#1a3692', \
              '#000066']

    plt.figure()
    plt.tight_layout()

    # if plot_dark_of_best:
    #     #make colors match
    #     color_cycle = plt.rcParams['axes.color_cycle']
    #     new_cycle = []
    #     for x in color_cycle:
    #         new_cycle.append(x)
    #         new_cycle.append(x)
    #     plt.rcParams['axes.color_cycle'] = new_cycle

    c = 0 # set index to 0 for color iteration
    
    for iv_file in champs['path']:
        iv_path = fix_path(iv_file)

        if iv_path:

            v, j = loadxIV(iv_path)
            log_file = champs[champs.path == iv_file]['logF'].values[0]
            # if log_file.find('corrected')>0:
            #     j*= 0.0919/0.1

            label_ = champs[champs.path == iv_file]['bin'].values[0]
            ls = '-'
            print(label_, iv_path)   # To find the best devices
            
            if not plot_dark_of_best and label_ == '$\mathregular{C_{60}}$':
                ls = '--'

            if plot_light_of_best:
                if logplot:
                    plt.semilogy(v, abs(j), label = label_, linewidth = 5, linestyle = ls , color = colors[c])
                else:
                    plt.plot(v, j, label = label_, linewidth = 5, linestyle = ls , color = colors[c])
                    plt.ylim([-4.3, 4.3])

            if plot_dark_of_best:
                r = '0.div'
                dark_address = ""
                l = len(iv_path)
                for i in range(len(iv_path)):
                    j = i-l
                    if j not in range (-6, -1):
                        dark_address += iv_path[j]
                    else:
                        dark_address += r[j+6]
#                dark_address = iv_path.replace('liv', 'div') # does not work if e.g. 1.div1 doesn't exist
                v, j = loadxIV(dark_address)
                if logplot:
                    if not plot_light_of_best:
                        plt.semilogy(v, abs(j), label = label_+' (Dark)', linewidth = 4, linestyle = '--' , color = colors[c])
                    else:
                        plt.semilogy(v, abs(j), linewidth = 4, linestyle = '--' , color = colors[c])

                else:
                    if not plot_light_of_best:
                        plt.plot(v, j, label = label_+' (Dark)', linewidth = 4, linestyle = '--' , color = colors[c])
                        plt.ylim([-1, 5])
                    else:
                        plt.plot(v, j, linewidth = 4, linestyle = '--' , color = colors[c])

            c += 1

    # CHANGE THIS TO ADJUST AXES LIMITS

#    plt.ylim([-10, 8])
#    plt.xlim([-1,1.5])
    plt.xlim([0,1.5])

    plt.ylim([-8, 8])
#    plt.ylim([-10, 8])
#    plt.ylim([-2, 5])
#    plt.ylim([-1, 5])
#    plt.xlim([0,1.5])

#    plt.plot([1.45,1.45],[-5,5], '-', color='xkcd:grey', linewidth=1) #plot a line from (0.01,-10) to (0.01, 10) for axis
#    plt.plot([0.01,0.01],[-5.9,10], '-', color='xkcd:grey', linewidth=1) #plot a line from (0.01,-10) to (0.01, 10) for axis
#    plt.plot([0.01,0.01],[-5.9,10], '-', color='xkcd:grey', linewidth=1) #plot a line from (0.01,-10) to (0.01, 10) for axis
#    plt.plot([0,10],[-6,-6], '-', color='xkcd:black', linewidth=4) #plot a line from (-10, 5.9) to (10, 5.9) for axis
    
    plt.axhline(0, color='xkcd:black', linewidth=2)
#    plt.axvline(0, color='xkcd:black', linewidth=4)
    plt.minorticks_on()
    plt.xlabel('Voltage / V', fontsize = 15,fontweight='medium')
    plt.ylabel('Current Density / mA/$\mathregular{cm^2}$', fontsize = 15, fontweight='medium')

    plt.grid(False)
    plt.rcParams['figure.facecolor']='xkcd:white'
    plt.rcParams['figure.edgecolor']='xkcd:white'
    plt.tick_params(labelsize=15, direction='in', axis='both', which='major', length=8, width=2)
    plt.tick_params(labelsize=15, direction='in', axis='both', which='minor', length=4, width=2)
#    plt.box(on=None)

#    plt.legend(loc='upper left',prop={'size':20})
#    plt.legend(loc='lower right',prop={'size':20})
    plt.legend(loc='best',prop={'size':15})

    plt.show()

def getUserSettings(no = 0):
    '''UI'''
    success = 0
    while not success:
        if no == 0:
            dir_ = input('Please enter the directory (full path) where the log files are: ') #IMPORTANT: This is just the path, such as ~/Desktop, NOT including the log file name!
        elif no == 1:
            dir_ = input('Please enter the directory (full path) where the stats files are: ')
        dir_ = os.path.expanduser(dir_)
        print(dir_)
        if os.path.exists(dir_): 
            success = 1
        else: print('Path %s does not exist, please enter a valid path'%dir_)
    success = 0
    while not success:
        if no == 0:
            binFile = input('Please enter the csv file in which the sample descriptions (bins) are encoded: ')
        elif no == 1:
            binFile = input('Please enter the csv file in which the sample descriptions (stats) are encoded: ')
        if not binFile.endswith('.csv'):
            print('The specified file is not a csv. Please specify a csv following the template')
        else:
            success = 1
    return binFile, dir_

if __name__ == '__main__':
    print('\n\tWelcome to IV BoxPlots\n'+'--'*20)
    home = os.getcwd() #save starting dir
    binsFile, dir_ = getUserSettings(no=0)
    cull = input('Cull duplicates? (If the same sample has been measured several times only the highest efficiency will be kept) 0/1\n>')

    try:
        os.chdir(dir_)
        setsDf = parse(binsFile) #load binsFile as data frame
        
        logFiles = pd.unique(setsDf['logFileAlias'].ravel()) #make list of log files from bins csv file
        
        for f in logFiles:
            logfile_df = parse(f)
            for i in range(len(logfile_df)):
                processRow(dict(zip(logStruct, logfile_df.iloc[i])), f)
        if cull: df = cullDuplicates(df) #incase several were measured

        subPos =np.arange(4)+1#subplot positions
        titles = ['Jsc', 'Voc', 'Eff', 'FF']
        units = ['$mA.cm^-$$^2$', '$V$', '$\%$', '']
        units_terminal = ['mA/cm^2', 'V', '%', '']
        #plt.figure('test') Display an empty window ?

        # print(df)
        prune = input('Determine statistics of pruned data? (0/1)\n>')
        
        if int(prune):
            pruneData(df.sort_values(['bin']), binsFile)

        if boxplots:

            for i, t in enumerate(titles):
                #unfortunately subplot not currently possible with pandas boxplot when using by
                boxplot(df.sort_values(['bin']), t)
                plt.title(t)
                plt.xticks(rotation='0') #change the orientation of the x axis tick
                plt.ylabel(units[i])
                plt.show()

        if plot_best_IVs:
            plot_best_cells(df)

    except Exception:
        import traceback
        print(sys.exc_info())
        for frame in traceback.extract_tb(sys.exc_info()[2]):
            fname,lineno,fn,text = frame
            print("Error in %s on line %d" % (fname, lineno))
    finally: 
        os.chdir(home)
