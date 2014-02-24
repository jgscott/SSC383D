'''
An example of Gaussian process sampling.

Author: Wesley Tansey
Date: 2/12/14
'''
from matplotlib import cm
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import scoreatpercentile
from functools import partial
import csv

class GaussianProcess(object):
    '''
    A Gaussian Process is an infinite-dimensional extension of the multivariate
    Gaussian. It takes a mean function over vectors and a covariance function
    over vectors that returns a covariance matrix.
    '''
    def __init__(self, mean, cov):
        self.mean = mean
        self.cov = cov

    def sample(self, x):
        ''' Given a vector of x values, returns a vector of y samples. '''
        mu = self.mean(x)
        sigma = self.cov(x, x)
        return np.random.multivariate_normal(mu, sigma)

    def predict(self, query, x, y, noise, percentile=95):
        '''
        Given a vector of query x values, returns an 2xN matrix with the first
        row being the predicted y values and the second being the size of the
        95% confidence band.
        '''
        # The common weighting term for the mean and covariance updates.
        w = self.cov(query, x).dot(np.linalg.inv(self.cov(x,x) + noise * np.eye(x.shape[0])))
        # The posterior predictive mean.
        mu = w.dot(y)

        # The posterior predictive covariance.
        sigma2 = self.cov(query, query) - w.dot(self.cov(x, query))
        
        # Sample from the posterior a bunch of times
        samples = np.random.multivariate_normal(mu, sigma2, 1000)

        # Approximate the 95% confidence intervals (symmetric for MVN)
        upper = 100. - (100. - percentile) / 2.
        upper_ci = scoreatpercentile(samples, upper)
        lower_ci = scoreatpercentile(samples, 100. - upper)

        return np.array([mu, upper_ci, lower_ci])

    def pseudo_marginal_log_likelihood(self, x, y, noise):
        ''' The log-likelihood of the data given the model, with the constants removed. '''
        k_y = self.cov(x,x) + noise * np.eye(x.shape[0])
        return -np.log(np.linalg.det(k_y)) - y.dot(np.linalg.inv(k_y)).dot(y)

def squared_exponential(bandwidth, tau1, tau2, x1, x2):
    ''' The squared exponential kernel function. '''
    # Restructure x1 as a matrix with row duplicates.
    # Then subtract each column of x2 from each row.
    x = np.tile(x1, (x2.shape[0],1)).T - x2

    # Kronecker delta
    kron = x == 0
    
    # Squared exponential
    return tau1 * np.exp(-0.5 / bandwidth**2 * x**2) + tau2 * kron

def load_utilities():
    ''' Loads the utilities data from the CSV file. '''
    with open('utilities.csv', 'rb') as f:
        reader = csv.reader(f)
        reader.next()
        x = []
        y = []
        for line in reader:
            x.append(float(line[0]))
            y.append(float(line[1]) / float(line[2]))
        return (np.array(x), np.array(y))

def sample_gp_prior():
    ''' Experiment to sample from a Gaussian process prior and plot the results. '''
    # Parameters of the experiment
    NUM_POINTS = 100
    NUM_TRIALS = 4
    KERNEL = 'squared_exponential'
    COLORS = ['red', 'blue', 'green', 'gold']

    # Hyperparameters of the squared exponential function
    BANDWIDTH = 0.5
    TAU1 = 1.
    TAU2 = 0.5

    # Create a partial function so we can try multiple data points easily
    sqexp = partial(squared_exponential, BANDWIDTH, TAU1, TAU2)

    # Create our Gaussian process with zero mean and squared exponential covariance
    gp = GaussianProcess(lambda x: np.zeros(x.shape[0]), sqexp)

    # Initialize the plot
    ax = plt.axes([.1,.1,.8,.7])

    # Sample some data points along the x-axis
    x = np.random.random(NUM_POINTS)

    # Sample from the GP multiple times to see multiple trials on the same plot
    for i in xrange(NUM_TRIALS):

        # Sample response values for the x samples
        y = gp.sample(x)

        # Plot this trial
        plt.scatter(x, y, label='Trial {0}'.format(i+1), color=COLORS[i])

    # Pretty up the plot
    plt.xlabel('x')
    plt.ylabel('y')
    plt.figtext(.40,.9, 'Gaussian Process Simulation Results', fontsize=18, ha='center')
    plt.figtext(.40,.85, '{0} points, {1} kernel, bandwidth={2} , tau1={3}, tau2={4}'.format(NUM_POINTS, KERNEL, BANDWIDTH, TAU1, TAU2), fontsize=10, ha='center')
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=12)
    plt.savefig('gp_samples.pdf')
    plt.clf()

def utilities_posterior():
    ''' Perform Gaussian process regression on the data in utilities.csv. '''
    # Parameters of the experiment
    KERNEL = 'squared_exponential'
    COLORS = ['red', 'blue', 'green', 'gold']

    # Hyperparameters of the squared exponential function
    BANDWIDTH = 4.
    TAU1 = 0.5
    TAU2 = 0.

    # Variance of the observation noise
    NOISE = 1.

    # Create a partial function so we can try multiple data points easily
    sqexp = partial(squared_exponential, BANDWIDTH, TAU1, TAU2)

    # Create our Gaussian process with zero mean and squared exponential covariance
    gp = GaussianProcess(lambda x: np.zeros(x.shape[0]), sqexp)

    # Initialize the plot
    ax = plt.axes([.1,.1,.8,.7])

    # Load the utilities data from the 'utilities.csv' file
    x,y = load_utilities()

    # The query x values
    query = np.linspace(0,100,1001)

    # The predicted y values for the query
    predicted = gp.predict(query, x, y - y.mean(), NOISE) + y.mean()

    # Plot the observed data points
    plt.scatter(x, y, label='Observed data', color=COLORS[0])
    plt.plot(query, predicted[0], label='Posterior', color=COLORS[1])
    plt.fill_between(query, predicted[1], predicted[2], facecolor=COLORS[1], alpha=0.2)

    # Pretty up the plot
    plt.xlim(0,100)
    plt.xlabel('Avg. Monthly Temperature (Fahrenheit)')
    plt.ylabel('Avg. Daily Utility Bill ($)')
    plt.figtext(.40,.9, 'Gaussian Process Regression on Utilities Data', fontsize=18, ha='center')
    plt.figtext(.40,.85, '{0} points, {1} kernel, bandwidth={2} , tau1={3}, tau2={4}, noise={5}'.format(len(x), KERNEL, BANDWIDTH, TAU1, TAU2, NOISE), fontsize=10, ha='center')
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=12)
    plt.savefig('utilities_posterior.pdf')
    plt.clf()

def utilities_heatmap():
    ''' Enumerate the tau1 and bandwidth values to find the minimum of the negative log-likelihood '''
    # Fixed parameters of the experiment
    KERNEL = 'squared_exponential'
    COLORS = ['red', 'blue', 'green', 'gold']
    TAU2 = 0.
    NOISE = 0.5235608 # empirical bayes variance measurement

    # Range of the parameters to enumerate
    MIN_BANDWIDTH = 40.
    MAX_BANDWIDTH = 70.
    MIN_TAU1 = 5.
    MAX_TAU1 = 35.
    SAMPLES = 200

    # Load the utilities data from the 'utilities.csv' file
    x,y = load_utilities()
    y = y - y.mean()

    # Calculate the heatmap values of the log-likelihood
    log_likelihood = np.zeros((SAMPLES,SAMPLES))
    bandwidths = np.linspace(MIN_BANDWIDTH, MAX_BANDWIDTH, SAMPLES)
    tau1s = np.linspace(MIN_TAU1, MAX_TAU1, SAMPLES)
    for i,BANDWIDTH in enumerate(bandwidths):
        for j,TAU1 in enumerate(tau1s):
            # Create a partial function so we can try multiple data points easily
            sqexp = partial(squared_exponential, BANDWIDTH, TAU1, TAU2)

            # Create our Gaussian process with zero mean and squared exponential covariance
            gp = GaussianProcess(lambda x: np.zeros(x.shape[0]), sqexp)

            # Calculate the pseudo marginal log-likelihood
            log_likelihood[SAMPLES-j-1,i] = gp.pseudo_marginal_log_likelihood(x, y, NOISE)

    # Plot the results
    ax = plt.axes([.1,.1,.8,.7])
    heatmap = plt.imshow(log_likelihood, cmap=cm.coolwarm, extent=[MIN_BANDWIDTH, MAX_BANDWIDTH, MIN_TAU1, MAX_TAU1])
    plt.colorbar(heatmap)
    plt.xlabel('bandwidth')
    plt.ylabel('tau1')
    plt.figtext(.50,.9, 'Hyperparameter Enumeration on Utilities Data', fontsize=18, ha='center')
    plt.figtext(.50,.85, '{0} points, {1} kernel, bandwidth=[{2},{3}] , tau1=[{4},{5}], tau2={6}, noise={7}'.format(len(x), KERNEL, MIN_BANDWIDTH, MAX_BANDWIDTH, MIN_TAU1, MAX_TAU1, TAU2, NOISE), fontsize=10, ha='center')
    plt.savefig('utilities_heatmap.pdf')
    plt.clf()


if __name__ == '__main__':
    utilities_posterior() # regression
    utilities_heatmap() # MAP parameter discovery
    