'''

'''

from load_file_by_header_columns import load_file_by_headers
import numpy as np
import matplotlib.pyplot as plt
from IPython.display import set_matplotlib_formats
from enum import Enum
import math

def setup_environment():
    '''
        Call the normal things I call to get the matplotlib settings I like.
    '''
    prop_cycle = plt.rcParams['axes.prop_cycle']
    colors = prop_cycle.by_key()['color']
    set_matplotlib_formats('jpeg')
    plt.rcParams['savefig.dpi'] = 300
    plt.rcParams["figure.dpi"] = 300
    plt.rcParams['figure.autolayout'] = False
    plt.rcParams['font.family'] = 'serif'
    plt.rcParams['lines.linewidth'] = 1.0
    plt.rcParams['lines.markersize'] = 5
    params = {
        'figure.figsize': (4, 3)
    }
    plt.rcParams.update(params)

    SMALL_SIZE = 8
    MEDIUM_SIZE = 9
    BIGGER_SIZE = 11
    plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
    plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
    plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
    plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
    plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title
    return colors

ZERO_CONCENTRATION_WEIGHT = 1e-6
# if the loq is not changing, we can stop boot-strapping
SAME_LOQ_REL_THRESHOLD = 0.001
np.random.seed(1234)
setup_environment()

class optimize_type(Enum):
    LOD = 0
    LOQ = 1

def get_default_config():
    return {
        'grid_size': 100, 
        'max_boot_iters': 100,
        'min_boot_iters_for_measure': 10,
        'min_same_loq_count_for_accept': 25,
        'cv_threshold': 0.2,
        'debug_plot': False,
        'optimize_type': optimize_type.LOQ,
        'minimum_num_transitions': 5,
    }

def add_concentration_to_raw_data(raw_data_items, rep_to_conc_data):
    '''
        Each raw data item has a key 'replicate'.  Here we match that key
        with items from another dict, the rep_to_conc_data, which has
        'replicate' and 'analyte_concentration' keys
        Any item from a replicate without a concentration label is removed
    '''
    # make a dict of the rep to conc data file
    # use float because if any concentrations are missing, all of the concentrations will be strings
    rep_to_conc_map = {}
    for item in rep_to_conc_data:
        # use float because if any concentrations are missing, all of the concentrations will be strings
        conc = 0.0 if item['analyte_concentration']=='False' else item['analyte_concentration'].astype(np.float)
        rep_to_conc_map[item['replicate']] = {'analyte_concentration': conc, 'sample_type': item['sample_type']}
    
    filtered_items = []
    for idx,item in enumerate(raw_data_items):
        rep = item['replicate_name']
        if not rep in rep_to_conc_map:
            raise Exception('Could not find replicate {0} for item[{1}] = {2}'.format(rep, idx, item))
        # don't include any solvent points
        if rep_to_conc_map[rep]['sample_type'] != 'Solvent':
            concentration = rep_to_conc_map[rep]['analyte_concentration']
            if concentration >= 0.0:
                item['analyte_concentration'] = rep_to_conc_map[rep]['analyte_concentration']
                filtered_items.append(item)
    
    if not filtered_items:
        raise Exception('There are no items with a concentration label')

    return filtered_items

def normalize_reps_by_global_standards(data):
    '''
        If there are any 'standard_type' == 'Global Standard', they must be present for each replicate,
        and we'll add them all up for each replicate, and divide the area for each item by that global standard sum
    '''
    global_standards_by_rep = {}
    found_global_standard = False
    for item in data:
        if item['standard_type'] == 'Global Standard':
            found_global_standard = True
            rep = item['replicate_name']
            if not rep in global_standards_by_rep:
                global_standards_by_rep[rep] = []
            global_standards_by_rep[rep].append(item['area'])

    if found_global_standard:
        for item in data:
            rep = item['replicate_name']
            if not rep in global_standards_by_rep:
                raise Exception('Replicate {0} was not found in the global standard replicates'.format(rep))
            global_mean_area = np.mean(global_standards_by_rep[rep])
            if global_mean_area > 0.0:
                item['area'] /= global_mean_area
    return data

def map_data_by_peptide(data):
    '''
        take the list of items, and organize it by peptide sequence
        return a dict with peptide+protein sequence keys, and list values of the raw data items
    '''
    peptide_to_data_map = {}
    for item in data:
        peptide_prot = item['peptide_modified_sequence'] + item['protein_name']
        if not peptide_prot in peptide_to_data_map:
            peptide_to_data_map[peptide_prot] = []
        peptide_to_data_map[peptide_prot].append(item)
    return peptide_to_data_map

def load_skyline_file_and_map_concentrations(transition_area_file, rep_to_conc_file):
    '''
        transition_area_file : a file with at least these headers, in any order:
            Peptide Modified Sequence, Fragment Ion, Area, Standard Type, Replicate Name
        rep_to_conc_file : a file with at least these headers :
            Replicate, Analyte Concentration, Sample Type
        return a list of dictionaries with the keys
            'peptide_modified_sequence', 'precursor_mz', 'precursor_charge', 'fragment_ion', 'product_mz', 'product_charge', 'area', 
            'standard_type', 'replicate_name', 'protein_name', 'analyte_concentration'
    '''
    # load the skyline file
    data = load_file_by_headers(
        transition_area_file, 
        # these headers should be present (can have different capitalization, spaces, and no '_')
        ['peptide_modified_sequence', 'precursor_mz', 'precursor_charge',
            'fragment_ion', 'product_mz', 'product_charge', 
            'area', 'standard_type', 'replicate_name', 'protein_name'],
        # functions to convert each column
        [str, float, int, str, float, int, float, str, str, str])
    
    rep_to_conc_data = load_file_by_headers(
        rep_to_conc_file, 
        ['replicate', 'analyte_concentration', 'sample_type'],
        [str, float, str]
    )

    # map the replicates to their concentrations, adding the 'analyte_concentration' key to each item
    data = add_concentration_to_raw_data(data, rep_to_conc_data)
    # if there are global standard replicates, normalize the areas by them
    data = normalize_reps_by_global_standards(data)
    # create a dictionary with the peptide sequence as the keys
    return map_data_by_peptide(data)

def sort_two_lists( list1, list2 ):
    '''
        sort by a first list, and sort the second list
        with reference to the first list
    '''
    list1, list2 = (list(x) for x in zip(*sorted(zip(list1, list2), key=lambda pair: pair[0])))
    return list1, list2

def organize_transition_data(peptide_data):
    '''
        Given a list of items
        dict {'peptide_modified_sequence', 'fragment_ion', 'area', 'standard_type', 'replicate_name', 'analyte_concentration'}
        create calibration data dictionary, that has the same keys as any input item except now the concentrations area a list,
        and the areas for each fragment are in an ndarray with rows for each concentration and columns for each fragment
    '''
    fragment_to_conc_area = {}
    for item in peptide_data:
        fragment_ion = '{0}_+{1}'.format(item['fragment_ion'], item['product_charge'])
        if not fragment_ion in fragment_to_conc_area:
            fragment_to_conc_area[fragment_ion] = {
                'product_mz': item['product_mz'], 
                'product_charge': item['product_charge'],
                'analyte_concentration': [], 'area': []}
        fragment_to_conc_area[fragment_ion]['analyte_concentration'].append(item['analyte_concentration'])
        fragment_to_conc_area[fragment_ion]['area'].append(item['area'])

    # get a list of the concentrations from the first fragment
    # and make sure the other transitions all have the same number of concentrations
    concentrations = None
    first_frag = None
    for frag,data in fragment_to_conc_area.items():
        frag_conc = data['analyte_concentration']
        if not concentrations:
            concentrations = sorted(frag_conc)
            first_frag = frag
            num_concentrations = len(concentrations)
        if num_concentrations != len(frag_conc):
            raise Exception('Num concentrations {0} in {1} dont match {2} in {3}'.format(len(frag_conc), frag, num_concentrations, first_frag))
    
    # shape the fragment data into an ndarray with concentrations on the rows and fragments across the columns
    fragment_areas = np.ndarray(shape=(num_concentrations, len(fragment_to_conc_area)))
    fragment_names, fragment_mz, fragment_z = [],[],[]
    for col_idx, (fragment_ion, frag_data) in enumerate(fragment_to_conc_area.items()):
        fragment_names.append(fragment_ion)
        fragment_mz.append(frag_data['product_mz'])
        fragment_z.append(frag_data['product_charge'])
        # sort the data for this fragment by concentration
        conc_for_frag, area_for_frag = sort_two_lists(frag_data['analyte_concentration'], frag_data['area'])
        if len(conc_for_frag) != num_concentrations:
            raise Exception('Number of concentrations in {0} does not match expected number {1} for {2}'.format(len(conc_for_frag), num_concentrations, item['peptide_modified_sequence']))
        # add to the ndarray
        for row_idx, area in enumerate(area_for_frag):
            fragment_areas[row_idx][col_idx] = area
    
    calibration_data = {
        'peptide_modified_sequence': item['peptide_modified_sequence'],
        'precursor_mz': item['precursor_mz'],
        'precursor_charge': item['precursor_charge'],
        'standard_type': item['standard_type'],
        'protein_name': item['protein_name'],
        'concentrations': concentrations,
        'areas': fragment_areas,
        'fragment_names': fragment_names,
        'product_mz': fragment_mz,
        'product_charge': fragment_z,
    }
    return calibration_data

def load_and_package_data(transition_area_file, rep_to_conc_file):
    '''
        Take the raw data generated for each peptide, and convert into a form that 
        is convenient for finding the loq of any transition or group of transitions,
        namely, where each peptide has a concentrations and areas lists, where the areas
        list is a 2d list of concentrations down the rows and fragments across the columns

        calibration_data = {
            'peptide_modified_sequence': item['peptide_modified_sequence'],
            'standard_type': item['standard_type'],
            'protein_name': item['protein_name'],
            'concentrations': concentrations,
            'fragment_names': fragment_names,
            'areas': fragment_areas,
        }
        transition_area_file : a file with at least these headers, in any order:
            Peptide Modified Sequence, Fragment Ion, Area, Standard Type, Replicate Name
        rep_to_conc_file : a file with at least these headers :
            Replicate, Analyte Concentration
    '''
    try:
        pep_prot_data = load_skyline_file_and_map_concentrations(transition_area_file, rep_to_conc_file)
        cal_data = []
        for pep_prot,peptide_data in pep_prot_data.items():
            cal_data.append(organize_transition_data(peptide_data))
    except Exception as e:
        e.args = ('Problem for {0}\n{1}'.format(pep_prot, e),)
        raise
    return cal_data

class long_statistics:
    '''
        This class contains methods for computing average and covariance on vary large data sets.
        The clever way to compute average and covariance for a large data set doesn't store all the data
        but just the data that is needed, which is the last value, and the count of how many value have
        been collected.  This can be seen for the case of computing the mean :
            y[n] = 1/N Sum[i=1,N] ( x[i] )
        Now factor out the last raw value, x[n]
            y[n] = 1/N Sum[i=1,N-1]( x[i] ) + 1/N x[n]
        The previous average is actually 1/N-1 Sum[i=1,N-1] ( x[i] )
        So we can say
            y[n] = (N-1) / N * 1/N Sum[i=1,N-1]( x[i] ) + 1/N x[n]
                 = (N-1) / N * y[n-1] + 1/N x[n]
    '''
    def __init__(self):
        self.count = 0
        self.x = 0.0
        self.x2 = 0.0
        
    def add(self, value):
        '''
            Add a new value
        '''
        self.count += 1
        if self.count == 1:
            self.x = value
            self.x2 = value**2
        else:
            inverse_n, old_factor = self.compute_long_coefficients(self.count)
            self.x = old_factor * self.x + inverse_n * value
            self.x2 = old_factor * self.x2 + inverse_n * value * value
            
    def compute_long_coefficients(self, count):
        '''
            Get the coefficients needed for computing the statistics
            These are two values that depend on the iteration number
        '''
        if count == 0:
            raise Exception('Count cannot be zero')
            
        inverse_n = 1.0 / count
        old_factor = (count - 1.0) * inverse_n
        
        return inverse_n, old_factor
    
    def get_count(self):
        return self.count
    
    def get_mean(self):
        return self.x
    
    def get_stats(self):
        return {'mean': self.get_mean(), 'std': self.get_std(), 'variance': self.get_variance()}
    
    def get_std(self):
        return np.sqrt(self.get_variance())
    
    def get_variance(self):
        return self.x2 - self.x**2
    
    def reset(self):
        self.count = 0
        self.x = 0.0
        self.x2 = 0.
        
def get_unique_conc_and_mean_areas(concentrations, areas):
    '''
        From equal-length lists of concentrations and areas,
        average the areas belonging to the same concentrations,
        and return a list of the unique concentrations and mean areas
    '''
    unique_conc = {}
    for conc,area in zip(concentrations, areas):
        if not conc in unique_conc:
            unique_conc[conc] = []
        unique_conc[conc].append(area)
    
    mean_areas = [np.mean(x) for _,x in unique_conc.items()]
    unique_conc = [conc for conc,_ in unique_conc.items()]
    # make sure they are sorted by concentration          
    return sort_two_lists(unique_conc, mean_areas)

def bilinear_eval(p, x):
    '''
        Evaluate a bilinear function for x with parameters p
        y = max(p[2], p[0]*x + p[1])
    '''
    return np.maximum(p[2], p[0]*x + p[1])

def bilinear_fit(x_offset, x_data, y_data):
    '''
        Simple way of doing the fit : 
        All data that are less than or equal to the x_offset get their areas
        averaged to determine the horizontal offset
        All data that are greater than the x_offset are fit to a line.
    '''
    x_subset,y_subset, weights_data,weights_sub = [],[],[],[]
    y_offset_data = []
    for x,y in zip(x_data, y_data):
        w = ZERO_CONCENTRATION_WEIGHT if x <= 0.0 else 1./x
        # using > ensures there is always at least one sample used for the offset
        if x > x_offset:
            x_subset.append(x)
            y_subset.append(y) 
            weights_sub.append(w)
        else:
            y_offset_data.append(y)
        weights_data.append(w)
            
    if len(x_subset) < 2:
        # the y offset describes the data set if there are too few points for fitting.  
        # Append the last y datum if it's there
        y_offset_data += y_subset
        p = [0., 0., np.mean(y_offset_data)]
    else:
        p = np.polyfit(x_subset, y_subset, 1, w=weights_sub)
        # the y offset is the mean of the data < x_offset
        p = np.append(p, np.mean(y_offset_data))
        
    # evaluate the error of the fit, over all of the data
    error = np.sum([((bilinear_eval(p, x) - y) * w)**2 for x,y,w in zip(x_data, y_data, weights_data)])

    return {'params': p,'error': error, 'baseline_std': np.std(y_offset_data)}
    
def fit_conc_vs_area(concentrations, areas):
    '''
        Fit a set of concentrations and areas to a bilinear fit
        return a dict of {
            'params': [a, b, c] parameters for maximum(c, a*x + b) equation
            'error': sum of squared errors from fitting
            'baseline_std': standard deviation of the data that were deemed to be in the baseline
        }
    '''
    # compute the means of each concentration
    unique_concentrations, mean_areas = get_unique_conc_and_mean_areas(concentrations, areas)
    # evaluate the bilinear fit using each possible concentration value as the baseline items
    fits = []
    for conc in unique_concentrations:
        # return the bilinear fit params, the error, and the baseline standard deviation
        fits.append(bilinear_fit(conc, unique_concentrations, mean_areas))

    # sort on the error, and then the best fit is the first one after sorting
    return sorted(fits, key=lambda x:x['error'])[0]

def compute_lod(fit, unique_concentrations):
    '''
        Calculate LOD using the bilinear fit 
        y = maximum(ax + b, c) and the baseline standard deviation
        For the area value that is c + baseline_std, find the corresponding
        concentration value
        return that concentration value
        The concentration returned will be clamped to the bounds of the non-zero
        concentrations that were measured
    '''
    (a,b,c), baseline_std = fit['params'],fit['baseline_std']
    # this is the area that corresponds to the loq
    lod_area = c + baseline_std
    largest_conc = max(unique_concentrations)
    # determine the lower bounds concentration for the lod
    smallest_nonzero_conc = largest_conc
    for conc in unique_concentrations:
        if conc > 0.0:
            smallest_nonzero_conc = min(conc, smallest_nonzero_conc)

    # sanity check on the slope
    if a==0.0:
        lod_conc = largest_conc
    else:
        lod_conc = (lod_area - b)/a
    # keep the loq within the bounds of the concentrations that were used
    lod_conc = max(smallest_nonzero_conc, min(lod_conc, largest_conc))
    
    return lod_conc

def compute_bootstrap_params(concentrations, areas):
    '''
        Bootstrapping is a technique for determining statistical confidence limits
        on a variable like mean and standard deviation.  
        https://en.wikipedia.org/wiki/Bootstrapping_(statistics)
        What we do is randomly sample from a set of already acquired data points,
        and can use these data to compute a set of best fit curves, giving
        a way to compute standard deviation and mean at values of concentration
        that weren't explicitly measured.
        Here we input a set of concentrations and areas, which are flattened out,
        that is, all replicates of the measurements are included.  In this function
        we simply randomly sample with replacment from these data and output num_boot_iters
        new data sets.
        return a set of parameters for fitting a bootstrapped set of data
    '''
    size_arrays = len(concentrations)
    new_conc, new_areas = np.ndarray(shape=(size_arrays)), np.ndarray(shape=(size_arrays))
    # get a set of random numbers from 0 to the size of the array, and do that until we have an array
    # the size of the original array
    random_indices = np.random.randint(0, size_arrays, size=(size_arrays))
    for row_idx in range(size_arrays):
        rand_index = random_indices[row_idx]
        new_conc[row_idx] = concentrations[rand_index]
        new_areas[row_idx] = areas[rand_index]
    # return the fit for the bootstrapped data set
    return fit_conc_vs_area(new_conc, new_areas)['params']
        
def compute_bootstrapped_loq(lod, concentrations, areas, max_concentration, config):
    '''
        lod : limit of detection
        concentrations : all concentrations corresponding to the areas replicates
        max_concentration : the largest concentration measured
    '''    
    # from lod to the highest concentration, search for a concentration that meets the cv requirement
    # create a set of evenly spaced concentrations from lod to highest concentration
    conc_grid = np.linspace(lod, max_concentration, config['grid_size'])
    # a list of objects that efficiently accumulate mean and standard deviation
    area_grid = [long_statistics() for _ in conc_grid]
    
    # use boot strapping to compute the loq
    num_iters_with_same_loq = 0
    last_loq,final_loq = max_concentration,max_concentration
    means,stds = None,None
    for idx in range(config['max_boot_iters']):
        # make a bootstrapped data set and return its bilinear fit
        p = compute_bootstrap_params(concentrations, areas)
        # update the grid areas with the bootstrapped fit
        for cidx,conc in enumerate(conc_grid):
            area_grid[cidx].add(bilinear_eval(p, conc))
        # after a minimum number of iterations, we can start computing the loq
        if idx > config['min_boot_iters_for_measure']:
            means = [x.get_mean() for x in area_grid]
            stds = [x.get_std() for x in area_grid]
            # default cv is a failing one, if the mean is 0
            cvs = [config['cv_threshold']*2 if y <= 0.0 else x/y for x,y in zip(stds,means)]
            # start at high concentration side and go down until cv threshold is reached
            # start with the loq at the maximum concentration as the fall back
            loq = max_concentration         
            for cv_idx in range(len(cvs)-1 , -1, -1):
                if cvs[cv_idx] > config['cv_threshold']:
                    break
                else:
                    # the loq is last one that met the condition
                    loq = conc_grid[cv_idx]
                    
            # check if the loq has changed
            if loq > 0.0 and np.abs(loq - last_loq) / loq < SAME_LOQ_REL_THRESHOLD:
                num_iters_with_same_loq += 1
            else:
                num_iters_with_same_loq = 0
                
            # save this loq for comparison with the last one
            last_loq = loq
            final_loq = loq
            # quit iterating if we've had to same loq for quite a while
            if num_iters_with_same_loq > config['min_same_loq_count_for_accept']:
                break
                
    return {
        'loq': final_loq, 'means': means, 'stds': stds, 'conc_grid': conc_grid,
    }

def compute_quantitative_limits(
    concentrations, 
    areas, 
    config):
    '''
        Given a list of concentrations and areas
        config = {
            grid_size : number of points in grid from lod to max concentration for figuring out loq,
            max_boot_iters : maximum number of iterations for bootstrapping
            min_boot_iters_for_measure : minimum number of iterations before starting to measure the loq
            min_same_loq_count_for_accept : if the loq stays the same for this many iterations, then stop bootstrapping early
            cv_threshold : loq is defined as the largest value that still satisfies std/mean < cv_threshold
        }
    '''
    best_fit = fit_conc_vs_area(concentrations, areas)
    # compute the lod
    unique_concentrations, mean_areas = get_unique_conc_and_mean_areas(concentrations, areas)
    lod = compute_lod(best_fit, unique_concentrations)
    loq_result = compute_bootstrapped_loq(lod, concentrations, areas, max(unique_concentrations), config)   
        
    if config['debug_plot']:
        fig,ax = plt.subplots(figsize=(4,3))
        
        ax.plot(concentrations, areas, marker='o', linewidth=0)
        ax.plot(unique_concentrations, mean_areas, marker='x', linewidth=0)
        
        ax.plot(unique_concentrations, [bilinear_eval(best_fit['params'], x) for x in unique_concentrations])
        ax.axhline(y=best_fit['params'][2], linestyle=':', color='brown')
        ax.axvline(x=lod, linestyle='--', color='magenta')
        ax.axvline(x=loq_result['loq'], linestyle=':', color='cyan')
        
        ax.set_xlabel('Concentration')
        ax.set_ylabel('Area')
        
        min_bounds = [x-y for x,y in zip(loq_result['means'], loq_result['stds'])]
        max_bounds = [x+y for x,y in zip(loq_result['means'], loq_result['stds'])]
        ax.fill_between(loq_result['conc_grid'], min_bounds, max_bounds, alpha=0.25, color='orange')
        ax.set_xscale('log')
        ax.set_yscale('log')        
    return {
        optimize_type.LOQ: loq_result['loq'], optimize_type.LOD: lod
    }
        
def get_other_opt_type(opt_type):
    if opt_type == optimize_type.LOQ:
        return optimize_type.LOD
    else:
        return optimize_type.LOQ
    
def optimize_transitions(peptide_info, config):
    '''
        The peptide info includes a 2d array of areas, with concentration down the rows,
        and fragment across the columns.  Here we compute the best set of transitions as follows:
        The quant limit (LOD or LOQ) is computed for each transition separately, and the top N
        are retained.  For the rest of the transitions, we do a test to include or not include them,
        depending on whether the quant limit is improved or not.
    '''
    concentrations,areas = peptide_info['concentrations'],peptide_info['areas']
    opt_type = config['optimize_type']
    other_opt_type = get_other_opt_type(opt_type)
    quant_limit_and_fragment_index = []
    num_concentrations, num_fragments = len(concentrations), len(areas[0])
    # keep track of the lowest limits for both types
    max_concentration = max(concentrations)
    lowest_limit = {
        optimize_type.LOQ: max_concentration, 
        optimize_type.LOD: max_concentration,
    }
    # check the limits for each of the fragments, by themselves
    for fidx in range(num_fragments):
        quant_limit = compute_quantitative_limits(concentrations, areas[:, fidx], config)
        quant_limit_and_fragment_index.append({
            'quant_limit': quant_limit, 'fragment_index': fidx,
        })
        # keep track of the lowest limits
        lowest_limit[optimize_type.LOD] = min(
            lowest_limit[optimize_type.LOD], quant_limit[optimize_type.LOD])
        lowest_limit[optimize_type.LOQ] = min(
            lowest_limit[optimize_type.LOQ], quant_limit[optimize_type.LOQ])        
    
    # if all of the limits of the chosen type are the max concentration, and the other one
    # isn't at the max, then start out using the best one of the other type
    if lowest_limit[opt_type] == max_concentration and lowest_limit[other_opt_type] < max_concentration:
        # sorting on the other quant limit
        quant_limit_and_fragment_index.sort(key=lambda x:x['quant_limit'][other_opt_type])
    else:
        # sort on the quant limit, low to high
        quant_limit_and_fragment_index.sort(key=lambda x:x['quant_limit'][opt_type])
    
    # start with the first ones of them that have less than the maximum for the LOQ.  
    # this array will accumulate the area of transitions that have been accepted
    accepted_areas = np.zeros(shape=(num_concentrations))
    accepted_fragment_indices = []
    for idx,item in enumerate(quant_limit_and_fragment_index[:config['minimum_num_transitions']]):
        quant_limit,fragment_index = item['quant_limit'][opt_type], item['fragment_index']
        # if we have accepted at least one transition, and the next ones are at the maximum
        # then break here, and test the transitions in the next loop
        if idx > 0 and quant_limit >= max_concentration:
            break
        # add the information for this fragment
        accepted_fragment_indices.append(fragment_index)
        for ridx in range(num_concentrations):
            accepted_areas[ridx] += areas[ridx, fragment_index]
    
    # after the last loop, we have at least one fragment on the list, and up to minimum_num_transitions
    # this is the benchmark quantitative limit
    optimized_quant_limit = compute_quantitative_limits(concentrations, accepted_areas, config)
    
    # start where we left off in the above loop
    start_index = min(len(accepted_fragment_indices), config['minimum_num_transitions'])
    
    # now for the rest of the transitions, test each one at a time with the already accepted one(s), and possibly add each one.
    rejected_items = []
    lowest_limit = {
        optimize_type.LOQ: max_concentration, 
        optimize_type.LOD: max_concentration,
    }
    for item in quant_limit_and_fragment_index[start_index:]:
        possible_new_areas = np.copy(accepted_areas)
        fragment_index = item['fragment_index']
        for ridx in range(num_concentrations):
            possible_new_areas[ridx] += areas[ridx, fragment_index]
        prospective_quant_limit = compute_quantitative_limits(concentrations, possible_new_areas, config)
        
        # accept this transition if it helped the result
        if prospective_quant_limit[opt_type] < optimized_quant_limit[opt_type]:
            optimized_quant_limit = prospective_quant_limit
            accepted_areas = possible_new_areas
            accepted_fragment_indices.append(fragment_index)
        else:
            # save the limits in case we don't have enough limits at the end of this
            rejected_items.append({'quant_limit': prospective_quant_limit, 'fragment_index': fragment_index})
            # keep track of the lowest limits
            lowest_limit[optimize_type.LOD] = min(
                lowest_limit[optimize_type.LOD], prospective_quant_limit[optimize_type.LOD])
            lowest_limit[optimize_type.LOQ] = min(
                lowest_limit[optimize_type.LOQ], prospective_quant_limit[optimize_type.LOQ]) 
            
    # if we still don't have enough transitions, for the case where there were transitions at the maximum limit
    if len(accepted_fragment_indices) < config['minimum_num_transitions'] and len(rejected_items) > 0:
        # again sort by the other type if all the rest of them are the max
        if lowest_limit[opt_type] == max_concentration and lowest_limit[other_opt_type] < max_concentration:
            # sort on the other type
            rejected_items.sort(key=lambda x:x['quant_limit'][other_opt_type])
        else:
            # sort the limits to add the best ones
            rejected_items.sort(key=lambda x:x['quant_limit'][opt_type])
        num_transitions_needed = config['minimum_num_transitions'] - len(accepted_fragment_indices)
        
        # add the needed transitions
        for item in rejected_items[:num_transitions_needed]:
            accepted_fragment_indices.append(item['fragment_index'])
            for ridx in range(num_concentrations):
                accepted_areas[ridx] += areas[ridx, item['fragment_index']]
        # calculate the final limit
        optimized_quant_limit = compute_quantitative_limits(concentrations, accepted_areas, config)
    return {
        'quant_limit': optimized_quant_limit,
        'accepted_fragment_indices': accepted_fragment_indices,
        'num_transitions': len(accepted_fragment_indices),
        'peptide_modified_sequence': peptide_info['peptide_modified_sequence'],
    }

def compute_quant_limit_for_all_transitions(peptide_info, config):
    '''
        Don't optimize any transitions, just compute the the quantitative limits
    '''
    concentrations,areas = peptide_info['concentrations'],peptide_info['areas']
    num_concentrations, num_fragments = len(concentrations), len(areas[0])
    
    summed_areas = np.zeros(shape=(num_concentrations))
    for fidx in range(num_fragments):
        for ridx in range(num_concentrations):
            summed_areas[ridx] += areas[ridx, fidx]
            
    return {
        'quant_limit': compute_quantitative_limits(concentrations, summed_areas, config),
        'num_transitions': num_fragments,
    }

def compute_quant_limit_for_selected_transitions(peptide_info, selected_transition_indices, config):
    '''
        Don't optimize any transitions, just compute the the quantitative limits for a selected set of indices
    '''
    concentrations,areas = peptide_info['concentrations'],peptide_info['areas']
    num_concentrations, num_fragments = len(concentrations), len(areas[0])
    
    summed_areas = np.zeros(shape=(num_concentrations))
    for fidx in selected_transition_indices:
        for ridx in range(num_concentrations):
            summed_areas[ridx] += areas[ridx, fidx]
            
    return {
        'quant_limit': compute_quantitative_limits(concentrations, summed_areas, config),
        'num_transitions': len(selected_transition_indices),
    }

def optimize_transitions_for_peptides(calibration_data_for_each_peptide, config):
    '''
        Call optimize_transitions on each of the items in load_and_package_data
    '''
    optimized_information = []
    num_peptides = len(calibration_data_for_each_peptide)
    report_stride = max(1, int(0.1 * num_peptides))
    for cidx, cal_data in enumerate(calibration_data_for_each_peptide):
        if cidx % report_stride == 0:
            print('[{0}/{1}] : {2}'.format(cidx+1, num_peptides, cal_data['peptide_modified_sequence']))
        optimized_information.append(optimize_transitions(cal_data, config))
    return optimized_information

def calculate_quant_limit_for_peptides(calibration_data_for_each_peptide, config):
    '''
        Calculate the loq and lod for all the transitions present for each peptide
    '''
    nonopt_information = []
    num_peptides = len(calibration_data_for_each_peptide)
    report_stride = max(1, int(0.1 * num_peptides))
    for cidx, cal_data in enumerate(calibration_data_for_each_peptide):
        if cidx % report_stride == 0:
            print('[{0}/{1}] : {2}'.format(cidx+1, num_peptides, cal_data['peptide_modified_sequence']))
        nonopt_information.append(compute_quant_limit_for_all_transitions(cal_data, config))
    return nonopt_information

def write_quant_limit_information(file_name, calibration_data_for_each_peptide, quant_limit_info):
    '''
        Write the peptide and loq/lod information to a file
    '''
    with open(file_name, 'w') as f:
        f.write('Peptide Modified Sequence\tProtein Name\tPrecursor Mz\tPrecursor Charge\tLOD\tLOQ\tNum Transitions\n')
        for pep,quant in zip(calibration_data_for_each_peptide, quant_limit_info):
            if pep['standard_type'] != 'Global Standard':
                f.write('{0}\t{1}\t{2:.4f}\t{3}\t{4:.4e}\t{5:.4e}\t{6}\n'.format(
                    pep['peptide_modified_sequence'],
                    pep['protein_name'],
                    pep['precursor_mz'],
                    pep['precursor_charge'],
                    quant['quant_limit'][optimize_type.LOD],
                    quant['quant_limit'][optimize_type.LOQ],
                    quant['num_transitions'],
                ))
            
def load_quant_limit_information(file_name):
    '''
        Load a file saved by write_quant_limit_information
    '''
    return load_file_by_headers(
        file_name, 
        # these headers should be present (can have different capitalization, spaces, and no '_')
        ['peptide_modified_sequence', 'protein_name', 'precursor_mz', 'precursor_charge', 'lod', 'loq', 'num_transitions'],
        # functions to convert each column
        [str, str, float, int, float, float, int])

def geom_spaced_series(minimum, maximum, steps):
    '''
        Create a geometrically spaced series
    '''
    const = math.exp(math.log(maximum / minimum) / steps)
    return [minimum * const**i for i in range(steps)]

def get_item_name(item):
    return item['peptide_modified_sequence'] + item['protein_name']

def make_pep_dict(info):
    pep_dict = {}
    for idx,item in enumerate(info):
        name = get_item_name(item)
        pep_dict[name] = {'item': item, 'index': idx}
    return pep_dict

def compare_quantitative_limits(info1, info2):
    '''
        Compare the quantitative limits, being careful to match the peptides
        to their names, in case the lists are not in the same order
    '''
    lod_ratios, loq_ratios, list1_indices, list2_indices = [],[],[],[]
    names = []
    dict1, dict2 = make_pep_dict(info1), make_pep_dict(info2)
    for idx1, (name, item1) in enumerate(dict1.items()):
        if not name in dict2:
            raise Exception('Could not find name {} in second dict'.format(name1))
        item2 = dict2[name]
        lod_ratios.append(item1['item']['lod'] / item2['item']['lod'])
        loq_ratios.append(item1['item']['loq'] / item2['item']['loq'])
        list1_indices.append(idx1)
        list2_indices.append(item2['index'])
        names.append(name)
    return {
        'lod': lod_ratios, 'loq': loq_ratios, 'index1': list1_indices, 'index2': list2_indices, 'names': names,
    }