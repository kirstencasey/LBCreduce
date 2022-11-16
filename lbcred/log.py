import sys
import logging
from astropy.logger import AstropyLogger
from astropy.utils.console import color_print
from tqdm import tqdm


class LBCreduceLogger(AstropyLogger):

    def _set_defaults(self, level='DEBUG'):
        """
        Reset logger to its initial state
        """
        # Remove all previous handlers
        for handler in self.handlers[:]:
            self.removeHandler(handler)

        # Set levels
        self.setLevel(level)

        # Set up the stdout handlers
        self.sh = StreamHandler()
        self.th = TqdmHandler()
        self.addHandler(self.sh)
        self.propagate = False

    def start_tqdm(self):
        self.removeHandler(self.sh)
        self.addHandler(self.th)

    def end_tqdm(self):
        self.removeHandler(self.th)
        self.addHandler(self.sh)


class StreamHandler(logging.StreamHandler):
    """
    A specialized StreamHandler that logs INFO and DEBUG messages to
    stdout, and all other messages to stderr.  Also provides coloring
    of the output, if enabled in the parent logger.
    """

    def emit(self, record):
        '''
        The formatter for stderr
        '''
        if record.levelno <= logging.INFO:
            stream = sys.stdout
        else:
            stream = sys.stderr

        if record.levelno < logging.INFO:
            color_print(record.levelname, 'magenta', end='', file=stream)
        elif record.levelno < logging.WARN:
            color_print(record.levelname, 'green', end='', file=stream)
        elif record.levelno < logging.ERROR:
            color_print(record.levelname, 'brown', end='', file=stream)
        else:
            color_print(record.levelname, 'red', end='', file=stream)
        record.message = "{0}".format(record.msg)
        print(": " + record.message, file=stream)


class TqdmHandler(logging.StreamHandler):

    def emit(self, record):
        if record.levelno <= logging.INFO:
            stream = sys.stdout
        else:
            stream = sys.stderr

        if record.levelno < logging.INFO:
            color_print(record.levelname, 'magenta', end='', file=stream)
        elif record.levelno < logging.WARN:
            color_print(record.levelname, 'green', end='', file=stream)
        elif record.levelno < logging.ERROR:
            color_print(record.levelname, 'brown', end='', file=stream)
        else:
            color_print(record.levelname, 'red', end='', file=stream)
        record.message = '{0}'.format(record.msg)
        tqdm.write(': ' + record.message)


logging.setLoggerClass(LBCreduceLogger)
logger = logging.getLogger('LBCreduceLogger')
logger._set_defaults()
