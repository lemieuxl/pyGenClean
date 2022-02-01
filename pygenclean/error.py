"""Exceptions raised by pyGenClean."""


class ProgramError(Exception):
    """An :py:class:`Exception` raised in case of a problem in the pipeline.

    Args:
        msg (str): the message to print to the user before exiting.

    """
    def __init__(self, msg: str = ""):
        self.message = str(msg)
        super().__init__()

    def __str__(self) -> str:
        return self.message
