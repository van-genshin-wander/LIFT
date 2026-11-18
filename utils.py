import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import torch
import torch.nn as nn
import ot as pot
import math
from tqdm import tqdm
import seaborn as sns
import pandas as pd


def wasserstein(
    x0: torch.Tensor,
    x1: torch.Tensor,
    reg: float = 0.05,
    power: int = 1,
    **kwargs,
) -> float:

    ot_fn = pot.emd2
    a, b = pot.unif(x0.shape[0]), pot.unif(x1.shape[0])
    if x0.dim() > 2:
        x0 = x0.reshape(x0.shape[0], -1)
    if x1.dim() > 2:
        x1 = x1.reshape(x1.shape[0], -1)
    M = torch.cdist(x0, x1)
    if power == 2:
        M = M**2
    ret = ot_fn(a, b, M.detach().cpu().numpy(), numItermax=1e7)
    if power == 2:
        ret = math.sqrt(ret)
    return ret

def wasserstein_with_weights(x0: torch.Tensor, m0, x1: torch.Tensor, m1, method="exact", reg=0.05, power=1) -> float:

    ot_fn = pot.emd2
    a, b = m0/m0.sum(), m1/m1.sum()

    x0 = x0.reshape(x0.shape[0], -1)
    x1 = x1.reshape(x1.shape[0], -1)
    M = torch.cdist(x0, x1, p=2)
    if power == 2:
        M = M ** 2
    dist = ot_fn(a, b, M.cpu().numpy())
    return math.sqrt(dist) if power == 2 else dist

def compute_uot_plans(X, t_train, delta=1, use_mini_batch_uot=False, group_number=5, draw=False):
    '''
    Compute the UOT plan and semi-coupling

    Parameters
    -------------------------------------
    X: list of 2D arrays
        the data at different time points
    t_train: list
        the time points
    delta: float
        growth regularization parameter
    use_mini_batch_uot: bool
        whether to use mini-batch UOT computation
    group_number: int
        the number of groups for mini-batch UOT computation
    draw: bool
        whether to draw the predicted and true marginals

    Returns
    -------------------------------------
    uot_plans: list of 2D arrays
        the UOT plans between consecutive time points
    gamma0_plans: list of 2D arrays
        the semi-coupling at source time points
    gamma1_plans: list of 2D arrays 
        the semi-coupling at target time points
    '''
    uot_plans = []
    gamma0_plans = []
    gamma1_plans = []


    total_action = 0
    for i in tqdm(range(len(t_train)-1), desc="Computing UOT plans..."):
        X_source = X[i]
        X_target = X[i+1]
        n_source, n_target = X_source.shape[0], X_target.shape[0] 
        norm_2_dist = pot.dist(X_source, X_target, metric='euclidean')

        cos_sq = np.cos(np.minimum(norm_2_dist / (2 * delta), np.pi/2))**2
        cost_matrix = -np.log(np.where(cos_sq == 0, 1e-10, cos_sq))


        if not use_mini_batch_uot:
            a = np.ones(n_source)
            b = np.ones(n_target)
            G = pot.unbalanced.mm_unbalanced(a, b, cost_matrix, reg_m=[1.0,1.0])
            total_action += 2*(delta**2)*pot.unbalanced.mm_unbalanced2(a, b, cost_matrix, reg_m=[1.0,1.0], returnCost='total')
        else:
            # Mini-batch UOT computation
            a = np.ones(n_source)
            b = np.ones(n_target)
            G = np.zeros((n_source, n_target))

            # Shuffle indices
            source_perm = np.arange(n_source)
            np.random.shuffle(source_perm)
            target_perm = np.arange(n_target)
            np.random.shuffle(target_perm)

            # Split indices into groups
            source_indices = np.array_split(source_perm, group_number)
            target_indices = np.array_split(target_perm, group_number)

            # for src_idx in source_indices:
            #     for tgt_idx in target_indices:
            for src_idx,tgt_idx in zip(source_indices,target_indices):
                    sub_cost_matrix = cost_matrix[np.ix_(src_idx, tgt_idx)]
                    sub_a = a[src_idx]
                    sub_b = b[tgt_idx]
                    G_sub = pot.unbalanced.mm_unbalanced(sub_a, sub_b, sub_cost_matrix, reg_m=[1.0,1.0])
                    G[np.ix_(src_idx, tgt_idx)] = G_sub


        gamma0_plan = ((a / G.sum(1))[:, None]) * G 
        gamma1_plan = (b/G.sum(0))*G

        # total_action += 2*(delta**2)*action_value(G, a, b, cost_matrix, reg_m=[1.0,1.0])

        uot_plans.append(G)
        gamma0_plans.append(gamma0_plan)
        gamma1_plans.append(gamma1_plan)

        print(((gamma1_plans[i]- gamma0_plans[i])<0).any())
        
        if draw:
            source_pred = G.sum(1)
            tar_pred = G.sum(0)
            
            fig = plt.figure(figsize=(15, 5))

            plt.subplot(131)
            plt.plot(a, label = f'source_true_{i}')
            plt.plot(source_pred, label = f'source_pred_{i}')
            plt.legend()

            plt.subplot(132)
            plt.plot(b, label = f'target_true_{i+1}')
            plt.plot(tar_pred, label = f'target_pred_{i+1}')
            plt.legend()
            
            plt.subplot(133)
            plt.scatter(X_source[:,0],X_source[:,1],s=source_pred*10, alpha=0.5)
            plt.show()

    print('true_action', total_action)

    return uot_plans, gamma0_plans, gamma1_plans, total_action


def sample_map(pi: np.ndarray, batch_size: int = 256, replace: bool = True):
    '''
    Sample (x0, x1) indices from the coupling pi(x0, x1)

    Parameters
    ---------------------------------
    pi: 2D array
        the coupling matrix
    batch_size: int
        the number of samples to draw
    replace: bool  
        whether to sample with replacement

    Returns
    ---------------------------------
    i_samples: 1D array (batch,)
        sampled indices from source
    j_samples: 1D array (batch,)
        sampled indices from target
    '''
    # Compute row sums and probabilities
    row_sums = pi.sum(axis=1)
    total_sum = row_sums.sum()
    row_probs = row_sums / total_sum  # Compute probabilities for i

    # Sample i according to row probabilities
    i_samples = np.random.choice(pi.shape[0], p=row_probs, size=batch_size, replace=replace)

    # Sample j for each sampled i
    j_samples = np.zeros(batch_size, dtype=int)
    for idx, i in enumerate(i_samples):
        # Compute sampling probabilities for j in the current row
        row_p = pi[i] / row_sums[i]  # Normalize to ensure probabilities sum to 1
        j_samples[idx] = np.random.choice(pi.shape[1], p=row_p)

    return i_samples, j_samples


def sample_from_ot_plan(ot_plan: np.ndarray, x0: torch.Tensor, x1: torch.Tensor, batch_size: int = 256):
    '''
    Sample (x0, x1) indices from the coupling pi(x0, x1)

    Parameters
    ---------------------------------
    ot_plan: 2D array (N, M)
        the coupling matrix
    x0: 2D array (N, dim)
        source samples
    x1: 2D array (M, dim)
        target samples
    batch_size: int
        the number of samples to draw

    Returns
    ---------------------------------
    x0[i]: 2D array (batch, dim)
        sampled source points
    x1[j]: 2D array (batch, dim)
        sampled target points
    i: 1D array (batch,)
        sampled source indices 
    j: 1D array (batch,)
        sampled target indices 
    '''
    i, j = sample_map(ot_plan, batch_size, replace=True)
    return x0[i], x1[j], i, j


