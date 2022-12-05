"""SUSI Logging module"""
import os
import logging


class Logging:
    """
    Facade for python logging.
    Provides easy setup and configuration methods.
    """
    #: Default name of the SUSI logger
    LOG_NAME = 'FDT'
    #: Default Log Format
    FORMAT = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')

    @staticmethod
    def init_console(level: int = logging.INFO) -> None:
        """
        Init console logging.
        """
        # Delete Jupyter notebook root logger handler, see https://github.com/ipython/ipython/issues/8282
        logger = logging.getLogger()
        logger.handlers = []

        # Create root logger with minimum supported log level
        logger = logging.getLogger(Logging.LOG_NAME)
        logger.setLevel(level)

        # create console logging handler
        ch = logging.StreamHandler()
        ch.setLevel(logging.DEBUG)
        ch.setFormatter(Logging.FORMAT)
        logger.addHandler(ch)

    @staticmethod
    def init_file(file_name: str, level: int = logging.INFO) -> None:
        """
        Init logging to a log file.
        """
        os.makedirs(os.path.dirname(file_name), exist_ok=True)
        logger = Logging.get_logger()
        fh = logging.FileHandler(file_name)
        fh.setLevel(level)
        fh.setFormatter(Logging.FORMAT)
        logger.addHandler(fh)

    @staticmethod
    def get_logger() -> logging.Logger:
        """
        Retrieve the SUSI default logger
        """
        return logging.getLogger(Logging.LOG_NAME)

    @staticmethod
    def set_log_level(level: int) -> None:
        """
        Set the logging level of the application. Per default the level is set to INFO (1).

        ### Params
         - level: (int) The level to set, where
            - 0: DEBUG
            - 1: INFO
            - 2: WARNING
            - 3: ERROR
            - 4: CRITICAL
        """
        log = Logging.get_logger()
        if level < 1:
            log.setLevel(logging.DEBUG)
        elif level < 2:
            log.setLevel(logging.INFO)
        elif level < 3:
            log.setLevel(logging.WARNING)
        elif level < 4:
            log.setLevel(logging.ERROR)
        else:
            log.setLevel(logging.CRITICAL)

    @staticmethod
    def welcome(args, script_name, version: str = None, log_dir=None) -> None:
        """
        Print a welcome message in CLI scripts.
        """
        log = Logging.get_logger()
        log.info("========== %s ==============", script_name)
        if version is not None:
            log.info('Version: %s', version)
        log.info('Arguments: %s', args)
        if log_dir is not None:
            log.info('Logging directory: %s', log_dir)
        log.info("======================================")
