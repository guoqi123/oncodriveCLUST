"""
Contains OncodriveCLUSTLError class
"""


class OncodriveCLUSTLError(Exception):
    """Base class for exceptions in OncodriveCLUSTL module"""
    pass


class UserInputError(OncodriveCLUSTLError):
    """
    Exception raised for errors in the input
    """
    def __init__(self, message):
        """
        Args:
            message (str): error found

        Returns:
            None
        """
        self.message = message
