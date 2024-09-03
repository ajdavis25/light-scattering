import logging


def setup_logging():
    """
    configure logging settings
    """
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
