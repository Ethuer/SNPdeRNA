import numpy as np

# classes/functions go here
def viterbi_alg(A_mat, O_mat, observations):
    # get number of states
    num_obs = observations.size
    num_states = A_mat.shape[0]
    # initialize path costs going into each state, start with 0
    log_probs = np.zeros(num_states)
    # initialize arrays to store best paths, 1 row for each ending state
    paths = np.zeros( (num_states, num_obs+1 ))
    paths[:, 0] = np.arange(num_states)
    # start looping
    for obs_ind, obs_val in enumerate(observations):
        # for each obs, need to check for best path into each state
        for state_ind in xrange(num_states):
            # given observation, check prob of each path
            temp_probs = log_probs + \
                         np.log(O_mat[state_ind, obs_val]) + \
                         np.log(A_mat[:, state_ind])
            # check for largest score
            best_temp_ind = np.argmax(temp_probs)
            # save the path with a higher prob and score
            paths[state_ind,:] = paths[best_temp_ind,:]
            paths[state_ind,(obs_ind+1)] = state_ind
            log_probs[state_ind] = temp_probs[best_temp_ind]
    # we now have a best stuff going into each path, find the best score
    best_path_ind = np.argmax(log_probs)
    # done, get out.
    return (best_path_ind, paths, log_probs)
 
# main script stuff goes here
if __name__ == '__main__':
    # the transition matrix
    A_mat = np.array([[.6, .4], [.2, .8]])
    # the observation matrix
    O_mat = np.array([[.5, .5], [.15, .85]])
    # sample heads or tails, 0 is heads, 1 is tails
    num_obs = 20
    observations1 = np.random.randn( num_obs )
    
    observations1[observations1>0] = 1
    observations1[observations1<=0] = 0
##    print observations1    <<<---   enter the initial SNP calling  Bernoulli probabilities
    # we have what we need, do viterbi
    best_path_ind, paths, log_probs = viterbi_alg(A_mat, O_mat, observations1)
##    print best_path_ind, paths, log_probs
    print "obs1 is " + str(observations1)
    print "obs1, best path is" + str(paths[best_path_ind,:])



    
##    # change observations to reflect messed up ratio
##    observations2 = np.random.random(num_obs)
##    observations2[observations2>.85] = 0
##    observations2[observations2<=.85] = 1
##    # majority of the time its tails, now what?
##    best_path_ind, paths, log_probs = viterbi_alg(A_mat, O_mat, observations2)
##    print "obs2 is " + str(observations1)
##    print "obs2, best path is" + str(paths[best_path_ind,:])
##    #have it switch partway?
##    best_path_ind, paths, log_probs = viterbi_alg(A_mat, \
##                                          O_mat, \
##                                          np.hstack( (observations1, observations2) ) )
##    print "obs12 is " + str(np.hstack( (observations1, observations2) ) )
##    print "obs12, best path is" + str(paths[best_path_ind,:])

