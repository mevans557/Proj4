import numpy as np
import scipy as sp


from graph_handling import draw_lattice


# For passing to draw_lattice() to easily convert the multiple state vectors
# into one structure for memory efficiency
state_dict = {0: 'H', 1: 'S', 2: 'I'}


def naive_iterateSI(M, Si, Ii):
    '''
    Naively time-step iterates through graph-based SI
    Params:
    M: An adjacency matrix of a graph
    Si: Initial susceptible population vector
    Ii: Initial infected population vector
    Returns a tuple of time, infected, susceptible arrays
    '''
    # Track populations over time
    Is = np.array([np.sum(Ii)])
    Ss = np.array([np.sum(Si)])
    ts = np.array([0])
    while 1 in Si:
        # Ensure exit after no change
        oldIi = np.copy(Ii)
        # Set up state vectors to display correctly. NOTE: Only used for display
        statevec = Si + 2*Ii
        states = {x: state_dict[statevec[x]] for x in range(Ii.size)}
        # Infect neighbours
        Ii = Ii + ((M @ Ii) * Si)
        # These style of expression reset all entries to 1 or 0
        Ii[Ii > 1] = 1
        Ii[Ii < 0] = 0
        # Make sure vectors sum to 1
        Si = Si - Ii
        Si[Si < 0] = 0
        Si[Si > 1] = 1
        # Track populations
        ts = np.append(ts, ts[-1]+1)
        Is = np.append(Is, np.sum(Ii))
        Ss = np.append(Ss, np.sum(Si))
        # Draw the graph
        # draw_lattice(M, states, np.sqrt(Ii.size).astype(int))
        # break loop if no updates occurred
        if (np.array_equal(Ii, oldIi)):
            break
    statevec = Si + 2*Ii
    states = {x: state_dict[statevec[x]] for x in range(Ii.size)}
    # draw_lattice(M, states, np.sqrt(Ii.size).astype(int))
    # Return the tracked data
    return ts, Is, Ss


def Gillespie_iterateSI(M, Si, Ii, kI):
    '''
    Uses Gillespie-style algorithm to iterate through graph-based SI
    Params:
    M: An adjacency matrix of a graph
    Si: Initial susceptible population vector
    Ii: Initial infected population vector
    kI: Probability of an infection on an edge in some time unit
    Returns a tuple of time, infected, susceptible arrays
    '''
    # Draw initial graph
    statevec = Si + 2*Ii
    states = {x: state_dict[statevec[x]] for x in range(Ii.size)}
    # draw_lattice(M, states, np.sqrt(Ii.size).astype(int))
    # Track populations over time
    Is = np.array([np.sum(Ii)])
    Ss = np.array([np.sum(Si)])
    ts = np.array([0])
    while 1 in Si:
        # Ensure exit after no change
        oldIi = np.copy(Ii)
        # Set up state vectors to display correctly. NOTE: Only used for display
        statevec = Si + 2*Ii
        states = {x: state_dict[statevec[x]] for x in range(Ii.size)}
        # Get infectable nodes
        Iable = (M @ Ii) * Si
        # Calculate next reaction time
        if (np.sum(Iable) == 0):
            break
        alpha = kI * np.sum(Iable)
        r1 = np.random.rand()
        tau = 1/alpha * np.log(1/r1)
        # Select a reaction to occur
        change = np.random.choice(np.nonzero(Iable)[0])
        Ii[change] = 1
        Si[change] = 0
        # Track populations
        ts = np.append(ts, ts[-1] + tau)
        Is = np.append(Is, np.sum(Ii))
        Ss = np.append(Ss, np.sum(Si))
        # Draw the graph
        draw_lattice(M, states, np.sqrt(Ii.size).astype(int))
        # break loop if no updates occurred
        if (np.array_equal(Ii, oldIi)):
            break
    statevec = Si + 2*Ii
    states = {x: state_dict[statevec[x]] for x in range(Ii.size)}
    draw_lattice(M, states, np.sqrt(Ii.size).astype(int))
    # Return the tracked data
    return ts, Is, Ss


def Gillespie_iterateSIS(M, Si, Ii, kI, kS):
    '''
    Uses Gillespie-style algorithm to iterate through graph-based SI
    Params:
    M: An adjacency matrix of a graph
    Si: Initial susceptible population vector
    Ii: Initial infected population vector
    kI: Probability of an infection on an edge in some time unit
    kS: Probability of a node going back from I -> S in some time unit
    Returns a tuple of time, infected, susceptible arrays
    '''
    # Draw initial graph
    # statevec = Si + 2*Ii
    # states = {x: state_dict[statevec[x]] for x in range(Ii.size)}
    # draw_lattice(M, states, np.sqrt(Ii.size).as_type(int))
    # Track populations over time
    Is = np.array([np.sum(Ii)])
    Ss = np.array([np.sum(Si)])
    ts = np.array([0])
    while (ts[-1] < 100):
        if (np.count_nonzero(Ii) == 0):
            break
        # Set up state vectors to display correctly. NOTE: Only used for display
        # statevec = Si + 2*Ii
        # states = {x: state_dict[statevec[x]] for x in range(Ii.size)}
        # Get infectable nodes
        Iable = (M @ Ii) * Si
        # Calculate next reaction time
        alpha = (kI * np.sum(Iable)) + (kS * np.sum(Ii))
        r1 = np.random.rand()
        tau = 1/alpha * np.log(1/r1)
        # Select a reaction to occur
        r2 = np.random.rand()
        if (r2 < (1/alpha * kI*np.sum(Iable))):
            change = np.random.choice(np.nonzero(Iable)[0])
            Ii[change] = 1
            Si[change] = 0
        else:
            change = np.random.choice(np.nonzero(Ii)[0])
            Ii[change] = 0
            Si[change] = 1
        # Track populations
        ts = np.append(ts, ts[-1] + tau)
        Is = np.append(Is, np.sum(Ii))
        Ss = np.append(Ss, np.sum(Si))
        # Draw the graph
        # draw_lattice(M, states, np.sqrt(Ii.size).astype(int))
    statevec = Si + 2*Ii
    states = {x: state_dict[statevec[x]] for x in range(Ii.size)}
    draw_lattice(M, states, np.sqrt(Ii.size).astype(int))
    # Return the tracked data
    return ts, Is, Ss
