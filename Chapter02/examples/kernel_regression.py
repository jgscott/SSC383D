'''
An example of kernel regression on a noisy univariate function using a Gaussian
kernel.

Author: Wesley Tansey
Date: 2/5/14
'''
import matplotlib.pyplot as plt
import numpy as np

# The number of observed data points
NUM_POINTS = 50

# The number of testing data points
NUM_TEST = 10

# The parameters of our noisy observations
X_MIN = 0.
X_MAX = 10.
Y_MEAN = 5.
SIN_FREQ = 1.
NOISE_STDEV = 0.1

#BANDWIDTHS = np.array([0.01, 0.1, 0.25, 0.5, 1.0, 2.0])
BANDWIDTHS = np.array([0.1, 0.25, 0.5])
COLORS = ['green', 'orange', 'pink', 'purple', 'brown', 'yellow', 'black'] # max 8 functions

def regress(query, x, y, kernel, bandwidth):
    ''' Perform kernel regression and predict the value of y at query '''
    y_mean = y.mean()
    y_adj = y - y_mean
    w = weight(x, query, kernel, bandwidth)
    results = np.zeros(len(query))
    for i,q in enumerate(query):
        results[i] = np.sum(w[i] * y) / w[i].sum()
    return results

def weight(x1, x2, kernel, bandwidth):
    ''' x1 and x2 are vectors; len(x1) != len(x2) '''
    results = np.zeros((len(x2), len(x1)))
    for i,x in enumerate(x2):
        results[i] = 1. / bandwidth * kernel((x1-x)/bandwidth)
    return results

def gaussian_kernel(x):
    return 1. / np.sqrt(2. * np.pi) * np.exp(-x**2/2.)

def simulate_noisy_sin(x, mean, freq, noise_stdev):
    y = np.sin(x*freq)+mean
    if noise_stdev > 0:
        y += np.random.normal(0, noise_stdev, x.shape[0])
    return y

def mean_squared_error(x1, x2):
    return ((x1-x2)**2).mean()

def oos_validation(train, test, kernel, bandwidths):
    x = train[0]
    y = train[1]
    query = test[0]
    actual = test[1]
    predictions = np.array([regress(query, x, y, kernel, h) for h in bandwidths])
    mse = np.array([mean_squared_error(actual, prediction) for prediction in predictions])
    return (predictions, mse)

def sample_data(num_points):
    x = np.random.uniform(X_MIN, X_MAX, num_points)
    y = simulate_noisy_sin(x, Y_MEAN, SIN_FREQ, NOISE_STDEV)
    return (x,y)

def fitting_example():
    ''' An example of running a kernel regression to fit a dataset from some nonlinear function. '''
    # The parameters of our regression
    KERNEL = gaussian_kernel
    BANDWIDTH = 0.5

    # Simulate some noisy observations
    x,y = sample_data(NUM_POINTS)

    # The points to sample along the x-axis
    query = np.linspace(X_MIN, X_MAX, 501)

    # Get the regression estimate
    prediction = regress(query, x, y, KERNEL, BANDWIDTH)

    # Get the real curve
    actual = simulate_noisy_sin(query, Y_MEAN, SIN_FREQ, 0)

    # Plot the results
    ax = plt.axes([.1,.1,.8,.7])
    plt.scatter(x, y, label='Observations')
    plt.plot(query, actual, color='blue', label='Actual')
    plt.plot(query, prediction, color='red', label='Predicted')
    plt.xlabel('x')
    plt.ylabel('y')
    f_desc = 'sin(x*{0})+{1} + N(0,{2}^2)'.format(SIN_FREQ, Y_MEAN, NOISE_STDEV)
    plt.figtext(.40,.9, 'Linear Smoothing Results', fontsize=18, ha='center')
    plt.figtext(.40,.85, '{0} points, {1} kernel, {2} bandwidth, f(x)={3}'.format(NUM_POINTS, 'Gaussian', BANDWIDTH, f_desc), fontsize=10, ha='center')
    # Shink current axis by 20%
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    # Put a legend to the right of the current axis
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=12)
    plt.savefig('smooth.pdf')
    plt.clf()

def validation_example():
    ''' An example of measuring out-of-sample predictive validation for different bandwidths. '''
    # The parameters of our regression
    KERNEL = gaussian_kernel
    
    # Simulate some noisy observations
    train = sample_data(NUM_POINTS)

    # Simulate some noisy test data
    test = sample_data(NUM_TEST)

    # The points to sample along the x-axis to plot the function
    x_plot_points = np.linspace(X_MIN, X_MAX, 501)

    # The evaluated points using the estimators
    y_plot_points = np.array([regress(x_plot_points, train[0], train[1], KERNEL, h) for h in BANDWIDTHS])

    # Get the real curve
    actual = simulate_noisy_sin(x_plot_points, Y_MEAN, SIN_FREQ, 0)

    # Measure the MSE for each bandwidth choice
    mse = oos_validation(train, test, KERNEL, BANDWIDTHS)[1]

    # Plot the results
    ax = plt.axes([.1,.1,.8,.7])
    plt.scatter(train[0], train[1], label='Training data', color='blue')
    plt.scatter(test[0], test[1], label='Testing data', color='red')
    plt.plot(x_plot_points, actual, color='blue', label='Actual')
    for i,f in enumerate(y_plot_points):
        plt.plot(x_plot_points, y_plot_points[i], color=COLORS[i], label='h={0} (MSE: {1:.4f})'.format(BANDWIDTHS[i], mse[i]))
    plt.xlabel('x')
    plt.ylabel('y')
    f_desc = 'sin(x*{0})+{1} + N(0,{2}^2)'.format(SIN_FREQ, Y_MEAN, NOISE_STDEV)
    plt.figtext(.40, .9, 'Out-of-Sample Validation of Bandwidths', fontsize=18, ha='center')
    plt.figtext(.40, .85, '{0} points, {1} kernel, f(x)={2}'.format(NUM_POINTS, 'Gaussian', f_desc), fontsize=10, ha='center')
    # Shink current axis by 20%
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    # Put a legend to the right of the current axis
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=10)
    plt.savefig('bandwidths.pdf')
    plt.clf()



if __name__ == '__main__':
    fitting_example()
    validation_example()