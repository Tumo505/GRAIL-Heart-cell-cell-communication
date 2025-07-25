import numpy as np
from scipy import stats
from scipy.stats import rankdata

def van_elteren_test(data, group, strata):
    """
    Perform the Van Elteren test for stratified data.
    
    Parameters:
        data (array-like): The data values.
        group (array-like): Group labels (e.g., 0 or 1 for two groups).
        strata (array-like): Stratification labels (e.g., batch or chamber).
    
    Returns:
        float: Test statistic.
        float: p-value.
    """
    unique_strata = np.unique(strata)
    test_statistic = 0
    total_n = 0
    
    for stratum in unique_strata:
        # Subset data by stratum
        stratum_mask = strata == stratum
        stratum_data = data[stratum_mask]
        stratum_group = group[stratum_mask]
        
        # Rank data within the stratum
        ranks = rankdata(stratum_data)
        
        # Compute Wilcoxon rank-sum test statistic for the stratum
        group_1_ranks = ranks[stratum_group == 1]
        group_0_ranks = ranks[stratum_group == 0]
        n1 = len(group_1_ranks)
        n0 = len(group_0_ranks)
        total_n += n1 + n0
        
        # Sum ranks for group 1
        sum_ranks_group_1 = np.sum(group_1_ranks)
        
        # Compute test statistic for the stratum
        stratum_statistic = sum_ranks_group_1 - (n1 * (n1 + n0 + 1)) / 2
        test_statistic += stratum_statistic
    
    # Compute p-value (approximation using normal distribution)
    mean_statistic = total_n * (total_n + 1) / 4
    variance_statistic = total_n * (total_n + 1) * (2 * total_n + 1) / 24
    z_score = (test_statistic - mean_statistic) / np.sqrt(variance_statistic)
    p_value = 2 * (1 - stats.norm.cdf(abs(z_score)))
    
    return test_statistic, p_value