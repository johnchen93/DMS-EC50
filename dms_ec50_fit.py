import pandas as pd
import re
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from scipy.optimize import curve_fit
import sys
from helper import *
from plot import *

version = 3.1

def main():
    mergeGrowthPercentage()
    #plotSample()
    #plotSample()
    #plotSample()
    
    fitEC50(sample_step=None)
    
    #plotWorkingFit(sample_step=1,just_wt=True)
    # plotWorkingFit(sample_step=1,just_wt=False)
    # plotWorkingFit(sample_step=1,just_wt=False,warning='high end survival')
    
    cleanUpFits()
    
def mergeGrowthPercentage():
    # collect filenames
    data_dir = './data'
    f = []
    for (dirpath, dirnames, filenames) in os.walk(data_dir):
        f.extend(filenames)
        break
    naive = [x for x in f if not 'AMP' in x]
    sel = [ x for x in f if 'AMP' in x ]
    print(sel)
    
    # read the naive library to determine which variants truly exist
    # unlike usual the naive library is not used for scoring, only to see which variants are above sequencing noise before selection
    base_df = None
    filter_headers, freq_headers = [], []
    merge_headers = ['aa_mut','aa_mut_pos','identity','group']
    naive_conds = [x.split("_")[1] for x in naive]
    for i in range(len(naive)):
        file = naive[i]
        d = pd.read_csv(os.path.join(data_dir, file), sep='\t')
        d['aa_mut_pos'] = d.apply(lambda x: x['aa_mut_pos'] if x['aa_mut']!= 'WT' else f"G{x['group']}", axis=1)
        cond = naive_conds[i]
        header = f'>2x_noise_{cond}'
        freq_header = f'{cond}_freq'
        filter_headers.append(header)
        freq_headers.append(freq_header)
        if cond == 'rep1': # exclude the unselected groups with few reads
            d = d[ ~d['group'].isin([1,2]) ]
        d[header] = d.apply( lambda x: x['count']>5 and x['count']>2*x['expected_err_count'], axis=1)
        totals = d[d[header]].groupby('group')['count'].sum()
        print(totals)
        d[freq_header] = d.apply(lambda x: x['count']/totals[x['group']] if x[header] else np.nan, axis=1)
        d = d[merge_headers+[header, freq_header]]
        if base_df is None:
            base_df = d
        else:
            base_df = base_df.merge(d, on=merge_headers, how='outer')
    base_df[filter_headers] = base_df.fillna(value=False)[filter_headers]
    base_df['naive_reps_observed'] = base_df[filter_headers].sum(axis=1)
    base_df['naive_avg_freq'] = base_df[freq_headers].mean(axis=1)
    print(base_df['naive_avg_freq'].sum())
    real_vars = base_df[base_df['naive_reps_observed']>=2]
    print(real_vars)
    real_vars.to_csv('variants_above_noise.txt',sep='\t')
    
    # label groups that don't have enough reads to be excluded from the data
    # exclude G2-16-rep1, G2-128-rep1, G4-8-rep1, G6-8-rep1, G6-16-rep1, G6-32-rep1, G7-128-rep1 
    exclude = {'8.0_rep1':[4,6],'16.0_rep1':[2,6],'32.0_rep1':[6],'128.0_rep1':[2,7]}
    
    # read in the selected counts
    sel_headers = ['aa_mut','aa_mut_pos','group','identity'] # what to use to merge from the real variants
    d_list = []
    for file in sel:
        conc = file.split('_')[2]
        rep = file.split('_')[3][-1]
        
        d = pd.read_csv(os.path.join(data_dir, file),sep='\t', index_col=0)
        d['aa_mut_pos'] = d.apply(lambda x: x['aa_mut_pos'] if x['aa_mut']!= 'WT' else f"G{x['group']}", axis=1)
        print(file, conc, rep)
        exclude_groups = []
        for k, v in exclude.items():
            if k == f'{conc}_rep{rep}':
                exclude_groups = v
                print('excluding groups:', k, v,'from',file)
                ##print( d[(d['replicate']==int(rep))&(d['group'].isin(exclude_groups))] )
        
        
        # select only variants that exist in the naive library
        # if a variant exists in the naive library but is not in the selection data, add a dummy count of 1
        group_totals = d.groupby('group')['count'].sum().reset_index().rename(columns={'count':'total_count'}) # collect read counts before adding dummy scores
        d = d.merge(real_vars[sel_headers + ['naive_avg_freq']], on=sel_headers, how = 'right')
        d = d[~(d['group'].isin(exclude_groups))] # remove excluded groups
        print(d[d['group'].isin(exclude_groups)])
        d['dummy_count'] = d['count'].apply( lambda x: np.isnan(x) ) # mark dummy counts
        d['count'] = d['count'].fillna(1)
        d = d.merge(group_totals, on='group', how='left')
        
        # reapply some labels not available to missing variants
        d['replicate'] = int(rep)
        d['conc'] = float(conc)
        
        d_list.append(d)
    data = pd.concat(d_list)
    print(data[['group','replicate','conc']].drop_duplicates())

    # convert the counts to frequency
    data['freq'] = data['count']/data['total_count']
    #print(data)
    
    # merge the frequency with the growth after 6hr selection
    # read in the growth data
    d_list = []
    skip_step=2
    nrows = 8
    for i in range(1,3):
        d = pd.read_excel('VIM_library_OD600_post_6hr_sel.xlsx', sheet_name='OD growth', skiprows=i*skip_step+(i-1)*nrows, nrows=nrows-1 )
        d = d.melt(id_vars='group',var_name='conc',value_name='OD600_postsel')
        d['replicate'] = i
        d_list.append(d)
    growth = pd.concat(d_list)
    growth = growth[growth['conc']!='NO SEL']
    growth['conc'] = growth['conc'].astype(float)
    data=data.merge(growth, on=['group','replicate','conc'], how='left')
    
    # merge in the OD without selection for each group to estimate population size
    od = pd.read_excel('VIM_library_OD600_post_6hr_sel.xlsx', sheet_name='no sel OD')
    data = data.merge(od, on=['group','replicate'], how='left')
    
    # calculate the overall population size = freq in library * frac growth 
    data['pop_size'] = data['freq']* ( data['OD600_postsel'] * 10**9 )
    data['avg_naive_pop_size'] = data['naive_avg_freq'] * ( data['OD600_nosel'] * 10**9 )
    data['label'] = data.apply(lambda x: f'{str(x["aa_mut_pos"]).split(".")[0]}{x["aa_mut"]}', axis=1)
    
    data.to_csv('dms_survival.txt', sep='\t', index=False)

def plotSample():
    # plot a few sets from the data set to see what it looks like
    data = pd.read_csv('dms_survival.txt', sep='\t')
    data['label'] = data.apply(lambda x: f'{x["aa_mut_pos"]}{x["aa_mut"]}', axis=1)
    sel = data[['label']].drop_duplicates()
    print(sel.sample(50))
    
    sub = sel.sample(50)
    plot_data = data[data['label'].isin(sub['label'].values)]
    print(plot_data)
    order = plot_data.sort_values(['aa_mut_pos','aa_mut'])['label'].drop_duplicates()
    
    has_file = True
    i = 1
    outfile = ''
    while has_file:
        if os.path.isfile(f'pop_size_plot_sample_{i}.pdf'):
            i+= 1
        else:
            has_file = False
            outfile = f'pop_size_plot_sample_{i}.pdf'
    print(outfile)
    
    g = sns.FacetGrid(plot_data, col='label', hue='replicate', col_wrap=4, sharey=False, col_order = order)
    g.map( plt.scatter, 'conc', 'pop_size', facecolors='none' )
    print(g.axes.shape)
    for i in range(g.axes.shape[0]):
        ax = g.axes[i,]
        ax.set_ylim(bottom=0)
        ax.set_xlim(left=2)
        ax.set_xscale('log', basex=2)
    
    g.savefig(outfile, format='pdf', bbox_inches='tight')

# helper function for making initial guesses
def fitGuess(df):
    start_pop = df.at[df['conc'].idxmin(),'pop_size']
    end_pop = df.at[df['conc'].idxmax(),'pop_size']
    tot_diff = start_pop - end_pop
    warnings = []
    if tot_diff < 0: # ending population is high
        warnings.append( 'rising pop size' )
    if start_pop/end_pop < 10:
        warnings.append( 'low pop difference')
    df = df.reset_index(drop=True)
    ec50_guess = 0
    for i, row in df.iterrows():
        if i+1<len(df):
            next_row = df.loc[i+1,:].to_dict()
            diff = row['pop_size'] - next_row['pop_size']
            if diff > tot_diff/2:
                ec50_guess = (row['conc']+next_row['conc'])/2
                #ec50_guess = row['conc']
    if ec50_guess==0: # no guess made, no sharp drop
        ec50_guess = df['conc'].median()
    hill_guess = -1
    
    
    # check if the population is higher than the naive
    max_pop = df['avg_naive_pop_size'].values[0]
    top_max = max_pop
    top_min = start_pop/2
    if start_pop > max_pop:
        top_max = start_pop*1.2
    if start_pop < max_pop/2:
        ec50_guess = df['conc'].min()
        top_min = max_pop/2
    if start_pop < max_pop/5:
        warnings.append('immediate pop drop')
    if max_pop < 100000:
        warnings.append('low naive pop')
    if max_pop < 1000:
        warnings.append('very low naive pop')
    if end_pop > max_pop/4:
        warnings.append('high end survival')
    
    warnings = '|'.join(warnings)
    return {'top':start_pop, 'bot':end_pop, 'hill':hill_guess, 'ec50':ec50_guess,
            'topmax':top_max, 'botmax':np.inf, 'hillmax':np.inf, 'ec50max':512, 
            'topmin':top_min, 'botmin':end_pop*0.5, 'hillmin':-np.inf, 'ec50min':0, 
            'warnings':warnings}
            
def fitEC50(sample_step=100):
    # cycle through each variant, perform separate fits for each replicate
    # every 100 fits, plot one of the curves if successful
    # for every variant that fails, log the variant / replicate / and plot the curve so that a guess can be made
    data = pd.read_csv('dms_survival.txt', sep='\t')
    
    data = data.sort_values(['group','aa_mut_pos','aa_mut','replicate','conc'])
    print(data)
    guess_file = f'initial_guesses_v{version}.txt'
    fit_file = f'fitted_values_v{version}.txt'
    gb_headers = ['group','aa_mut_pos','aa_mut','replicate']
    
    if os.path.isfile(guess_file):
        guess_df = pd.read_csv(guess_file, sep='\t')
    else:
        guess_dict = {}
        for grouping, df in data.groupby(gb_headers):
            guess = fitGuess(df)
            guess_dict[grouping] = guess
        guess_df = pd.DataFrame.from_dict(guess_dict, orient='index').reset_index().rename(columns=dict(zip([f'level_{x}' for x in range(len(gb_headers))] ,gb_headers)))
        guess_df['rising_pop'] = guess_df.apply(lambda x: 'rising pop size' in x['warnings'], axis=1 )
        guess_df['low_pop_diff'] = guess_df.apply(lambda x: 'low pop difference' in x['warnings'], axis=1 )
        print(guess_df)
        guess_df.to_csv(guess_file,sep='\t', index=False)
    
    if os.path.isfile(fit_file):
        fit_df = pd.read_csv(fit_file, sep='\t')
    else:
        fit_dict = {}
        guess_headers = ['top','bot','hill','ec50']
        gdf = guess_df.groupby(gb_headers)
        for names, df in data.groupby(gb_headers):
            gu = gdf.get_group(names).reset_index(drop=True).loc[0,:].to_dict() # guess parameters
           
            p0 = [ gu[x] for x in guess_headers ] # guess parameter for input
            bounds = ( [gu[x+'min'] for x in guess_headers], [gu[x+'max'] for x in guess_headers] )
            vx, vy = df['conc'].values, df['pop_size'].values
            fitted = False
            try:
                popt, pcov = curve_fit(ECcurve, vx, vy, p0=p0, bounds=bounds) # curve fit
                fitted = True
            except:
                # no fit
                popt = [0 for x in range(4)]
                #print(f'Could not fit curve for {names}')
            fit_dict[names] = { 'top':popt[0],'bot':popt[1],'hill':popt[2],'ec50':popt[3], 'fitted':fitted, 'warnings':gu['warnings'] }
        fit_df = pd.DataFrame.from_dict(fit_dict, orient='index').reset_index().rename(columns=dict(zip([f'level_{x}' for x in range(len(gb_headers))] ,gb_headers)))
        fit_df.to_csv(fit_file, sep='\t', index=False)
    
    # plot all the fits that failed
    failed = fit_df[fit_df['fitted']==False]
    failed['warnings'] = failed.fillna('')['warnings']
    print(failed)
    
    if sample_step is None:
        return 
    outname = f'failed_fit_curves_v{version}.pdf'
    pp = PdfPages(outname)
    
    i = 0
    step_count = sample_step # plot every x failed plots
    pl_gb_headers = ['aa_mut_pos','aa_mut']
    pl = data.groupby(pl_gb_headers) # plotting data
    for names, df in failed.groupby(pl_gb_headers):
        i += 1
        if i%step_count!=0 and names[1]!='WT':
            continue
                
        plot_data = pl.get_group(names).sort_values(['replicate','conc'])
        plot_data = plot_data[plot_data['replicate'].isin(df['replicate'].unique())] # get just the replicate that failed
        
        g = sns.FacetGrid(plot_data, col='replicate', row='label', hue='dummy_count', palette = ['orange','grey'] )
        g.map( plt.scatter, 'conc', 'pop_size', facecolors='none' )
        g.map_dataframe(drawNaivePop)
        for row in range(g.axes.shape[0]):
            for col in range(g.axes.shape[1]):
                ax = g.axes[row,col]
                ax.set_ylim(bottom=0)
                ax.set_xlim(left=2)
                ax.set_xscale('log', basex=2)
                ax.text(1,0.95, '\n'.join(df['warnings'].values[0].split('|')), horizontalalignment='right', transform=ax.transAxes, fontsize='small', verticalalignment='top')
        g.savefig(pp, format='pdf', bbox_inches='tight')
        plt.close(plt.gcf())
    pp.close()

def plotWorkingFit(sample_step=50, just_wt=False, warning=''):
    fit_df = pd.read_csv(f'fitted_values_v{version}.txt', sep='\t')
    fit_df['warnings'] = fit_df.fillna('')['warnings']
    data = pd.read_csv('dms_survival.txt', sep='\t')
    data = data.sort_values(['group','aa_mut_pos','aa_mut','replicate','conc'])
    
    # plot the fits that worked, and the curve that was fitted
    worked = fit_df[fit_df['fitted']]
    if warning!='':
        worked = fit_df[fit_df['warnings'].str.contains(warning)]
    wtonly=''
    if just_wt:
        wtonly = '_wt_only'
        data = data[data['aa_mut']=='WT']
        worked = worked[worked['aa_mut']=='WT']
        sample_step = 1 # plot everything
    warn_label = ''
    if warning!='':
        warn_label='_'+warning
    outname = f'working_fit_curves_v{version}{wtonly}{warn_label}.pdf'
    pp = PdfPages(outname)
    pl_gb_headers = ['aa_mut_pos','aa_mut']
    pl = data.groupby(pl_gb_headers) # plotting data
    
    i = 0
    for names, df in worked.groupby(pl_gb_headers):
        i+=1
        if i%sample_step!=0:
            continue
        plot_data = pl.get_group(names).sort_values(['replicate','conc'])
        plot_data = plot_data[plot_data['replicate'].isin(df['replicate'].unique())] # get just the replicate that worked
        
        g = sns.FacetGrid(plot_data, col='replicate', row='label', hue='dummy_count', palette = ['orange','grey'] )
        g.map( plt.scatter, 'conc', 'pop_size', facecolors='none' )
        g.map_dataframe(sigmoidPlot, fits = df, xmin=2, xmax=256 )
        
        for row in range(g.axes.shape[0]):
            for col in range(g.axes.shape[1]):
                ax = g.axes[row,col]
                ax.set_ylim(bottom=0)
                ax.set_xlim(left=2)
                ax.set_xscale('log', basex=2)
        g.savefig(pp, format='pdf', bbox_inches='tight')
        plt.close(plt.gcf())
        i+=1
    pp.close()

def cleanUpFits():

    fits = pd.read_csv(f'fitted_values_v{version}.txt', sep='\t')
    fits['label']=fits.apply(lambda x: f'{str(x["aa_mut_pos"]).split(".")[0]}{x["aa_mut"]}_rep{x["replicate"]}', axis=1)
    working_qc = pd.read_excel(f'manual_fit_QC.xlsx', sheet_name='working fits' )
    working_qc['label']=working_qc.apply(lambda x: f'{str(x["aa_mut_pos"]).split(".")[0]}{x["aa_mut"]}_rep{x["replicate"]}', axis=1)
    failed_note = pd.read_excel(f'manual_fit_QC.xlsx', sheet_name='failed', skiprows=2 )
    failed_note['label']=failed_note.apply(lambda x: f'{str(x["aa_mut_pos"]).split(".")[0]}{x["aa_mut"]}_rep{x["replicate"]}', axis=1)
    print(failed_note)
    
    # ignore those with no points in the fitted transition, or those marked as notable form the failed set (don't know what to do with them yet)
    fits['ignore'] = (     ( fits['label'].isin(working_qc[~working_qc['missing_range'].isna()]['label']) ) 
                        |  ( fits['label'].isin(failed_note[failed_note['note']!='HIGH RESIST']['label']) )
                        |  ( fits['label'].isin(working_qc[~working_qc['exclude_other_note'].isna()]['label'] ) )   )
    # classify all failed variants with no other notes as dead, or those in the working fits who are also dead
    fits['dead'] = ( fits['label'].isin(working_qc[~working_qc['dead'].isna()]['label'] ))  |  ( (~fits['fitted'])  &  ( ~fits['label'].isin(failed_note['label']))  ) 
    # mark the high resistance ones
    fits['high_resist'] = ( fits['label'].isin(working_qc[~working_qc['high_resist'].isna()]['label']) )  |  ( fits['label'].isin(failed_note[failed_note['note']=='HIGH RESIST']['label']) )
    # fill in the variants with missing fits at the extremes of the experimental conditions
    fits['filled'] = fits.apply( lambda x: ( x['dead'] or x['high_resist'] ) and not x['ignore'], axis=1 )
    fits['filled_ec50'] = fits.apply( fillFit, axis=1 )
    
    print( working_qc[~working_qc['missing_range'].isna()] )
    print( fits[fits['high_resist']] )
    print( fits[fits['filled_ec50'].isna()]['ignore'].value_counts() )
    
    fits['log2_ec50'] = fits['filled_ec50'].apply( lambda x: np.log2(x) if not np.isnan(x) else np.nan )
    fits.to_csv('ec50_fitted_filled.txt', sep='\t', index=False)
    
def fillFit(row):
    min_ec, max_ec = 4, 256
    if row['dead']:
        return min_ec
    elif row['high_resist']:
        return max_ec
    else:
        return row['ec50']
        
def ECcurve(x, top, bottom, hill, ec50):
    y = bottom + (top-bottom)/(1+((ec50/x)**hill))
    return y

def logECcurve(log2x, top, bottom, hill, log2ec50):
    y = bottom + (top-bottom)/(1+2**((log2ec50-log2x)*hill) )
    return y

def sigmoidPlot(**kwargs):
    fits, xmin, xmax = kwargs['fits'], kwargs['xmin'], kwargs['xmax']
    if len(fits)==0:
        return
    ax = plt.gca()
    data = kwargs.pop("data")
    rep = data['replicate'].values[0]
    rep_df = fits[fits['replicate']==rep]
    if len(rep_df)==0:
        return
    popt = rep_df.iloc[0,:].to_dict()
    #print(popt)
    curve_x = np.linspace(xmin, xmax, 1000) # generate the curve for plotting
    curve_y = ECcurve(curve_x, popt['top'], popt['bot'], popt['hill'], popt['ec50'])
    ax.plot(curve_x, curve_y, zorder=-1, color='#ba953f')
    ax.axhline(data['avg_naive_pop_size'].values[0], color='red', linestyle='--', zorder=-1)

def drawNaivePop(**kwargs):
    ax = plt.gca()
    data = kwargs.pop("data")
    #print(data)
    ax.axhline(data['avg_naive_pop_size'].values[0], color='red', linestyle='--', zorder=-1)
    
def nextFileName(string):
    has_file = True
    i = 1
    outfile = ''
    while has_file:
        if os.path.isfile(string.format(i)):
            i+= 1
        else:
            has_file = False
            outfile = string.format(i)
        if i> 100: # something went wrong
            return None
    print(outfile)
    return outfile
    
if __name__ == '__main__':
    main()
