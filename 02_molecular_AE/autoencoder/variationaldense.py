from keras import backend as K
from keras.layers.core import Layer
from keras import initializations, activations

from seya.regularizers import GaussianKL


class VariationalDense(Layer):
    """VariationalDense
        Hidden layer for Variational Autoencoding Bayes method [1].
        This layer projects the input twice to calculate the mean and variance
        of a Gaussian distribution. During training, the output is sampled from
        that distribution as mean + random_noise * variance, during testing the
        output is the mean, i.e the expected value of the encoded distribution.
        Parameters:
        -----------
        batch_size: Both Keras backends need the batch_size to be defined before
            hand for sampling random numbers. Make sure your batch size is kept
            fixed during training. You can use any batch size for testing.
        regularizer_scale: By default the regularization is already properly
            scaled if you use binary or categorical crossentropy cost functions.
            In most cases this regularizers should be kept fixed at one.
    """
    def __init__(self, output_dim, batch_size, init='glorot_uniform',
                 activation='tanh',
                 weights=None, input_dim=None, regularizer_scale=1,
                 prior_mean=0, prior_logsigma=0, output_sample=False,
                 output_var=False,
                 **kwargs):
        self.prior_mean = prior_mean
        self.prior_logsigma = prior_logsigma
        self.regularizer_scale = K.variable(regularizer_scale)
        self.batch_size = batch_size
        self.init = initializations.get(init)
        self.activation = activations.get(activation)
        self.output_dim = output_dim
        self.initial_weights = weights
        self.input_dim = input_dim
        self.output_sample = output_sample
        self.output_var = output_var
        if self.input_dim:
            kwargs['input_shape'] = (self.input_dim,)
        self.input = K.placeholder(ndim=2)
        super(VariationalDense, self).__init__(**kwargs)

    def build(self):
        input_dim = self.input_shape[-1]

        self.W_mean = self.init((input_dim, self.output_dim))
        self.b_mean = K.zeros((self.output_dim,))
        self.W_logsigma = self.init((input_dim, self.output_dim))
        self.b_logsigma = K.zeros((self.output_dim,))

        self.trainable_weights = [self.W_mean, self.b_mean, self.W_logsigma,
                                  self.b_logsigma]

        self.regularizers = []
        reg = self.get_variational_regularization(self.get_input())
        self.regularizers.append(reg)

    def get_variational_regularization(self, X):
        mean = self.activation(K.dot(X, self.W_mean) + self.b_mean)
        logsigma = self.activation(K.dot(X, self.W_logsigma) + self.b_logsigma)
        return GaussianKL(mean, logsigma,
                          regularizer_scale=self.regularizer_scale,
                          prior_mean=self.prior_mean,
                          prior_logsigma=self.prior_logsigma)

    def get_mean_logsigma(self, X):
        mean = self.activation(K.dot(X, self.W_mean) + self.b_mean)
        logsigma = self.activation(K.dot(X, self.W_logsigma) + self.b_logsigma)
        return mean, logsigma

    def _get_output(self, X, train=False):
        mean, logsigma = self.get_mean_logsigma(X)
        if train or self.output_sample:
            # Temporary change, scale down size of noise
            if K._BACKEND == 'theano':
                eps = K.random_normal((X.shape[0], self.output_dim), std=self.regularizer_scale)
            else:
                eps = K.random_normal((self.batch_size, self.output_dim))
            # Temporary change, multiply by regularizer_scale
            return mean + self.regularizer_scale * K.exp(logsigma) * eps
        else:
            if self.output_var:
                return mean, logsigma
            else:
                return mean

    def get_output(self, train=False):
        X = self.get_input()
        return self._get_output(X, train)

    @property
    def output_shape(self):
        return (self.input_shape[0], self.output_dim)