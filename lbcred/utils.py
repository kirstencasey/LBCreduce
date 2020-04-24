import logging
from .log import logger


def setup_logger(level, log_fn=None):
    """
    Setup the pipeline logger.
    Parameters
    ----------
    level : str
        The log level (debug, info, warn, or error).
    log_fn : str (optional)
       Log file name.
    """
    if log_fn is not None:
        fh = logging.FileHandler(log_fn)
        formatter = logging.Formatter(
            '%(asctime)s | %(levelname)s: %(message)s',
            '%Y-%m-%d | %H:%M:%S')
        fh.setFormatter(formatter)
        logger.addHandler(fh)
    logger.setLevel(level.upper())

