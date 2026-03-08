import numpy as np
from matplotlib import pyplot as plt


def griddify(tin, tout, Iin):
    '''
    Params:
    tin: np array of non-grid times
    tout: np array of grid times
    Iin: np array of non-grid infecteds
    Returns a numpy array of grid infecteds, by sampling the non-grid Infecteds at grid times
    '''
    if (tin.size != Iin.size):
        print("mismatched size, zeros returned")
        return np.zeros(tout.size())

    Iout = []
    j = 0  # pointer to where we are in tin
    for t in tout:
        # If we are at the end, just append the last element
        if (j == tin.size-1):
            Iout.append(Iin[-1])
            continue

        # Not at the end, move j to the correct place and append
        while (t >= tin[j+1]):
            j += 1
            if (j >= tin.size-2):
                break
        Iout.append(Iin[j])

    return np.array(Iout)


if __name__ == "__main__":
    # test data
    testt = np.array([0, 0.12, 0.3, 0.5, 0.61, 0.73, 1])
    testI = np.array([10, 11, 12, 13, 12, 13, 14])
    grid = np.arange(0, 1.5, 0.2)

    # for testing purposes, plots original data correctly (step graph, not line)
    testt2 = testt - 0.001
    for i in range(testt2.size-1):
        testt2[i] = testt2[i+1]
    testt2[-1] = 1.1
    totaltestt = []
    totaltestI = []
    for i in range(testI.size):
        totaltestt.append(testt[i])
        totaltestt.append(testt2[i])
        totaltestI.append(testI[i])
        totaltestI.append(testI[i])

    gridI = griddify(testt, grid, testI)
    # plt.plot(testt, testI)
    plt.plot(totaltestt, totaltestI)
    plt.plot(grid, gridI)
    plt.show()
