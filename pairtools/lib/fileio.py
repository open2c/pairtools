import shutil
import subprocess
import shlex
import sys
import typing as tp
from dataclasses import dataclass

class ParseError(Exception):
    pass


# dictionary of allowed modes for each file type: 
# key: format, 
# value: list with allowed modes
PRESET_EXT2MODE = {
    'bam': ['w', 'r'],
    'gz': ['w', 'r', 'a'],
    'lz4': ['w', 'r', 'a'],
}

# dictionary of automatic opener commands
# key: tuple (file format, mode)
# value: dictionary with keys tools - 
#       'tool': name of the tool that the CommandRun will attempt to find via shutil,
#       'command': command line with the number of threads to be formatted before the execution.
# In summary:
# dict ('file_format', 'mode'): ['tool': 'tool_name', 'command': 'command formatted to set threads number']
PRESET_COMMANDS = {
    ('bam', 'w'): [
        {'tool': 'samtools', 'command': 'samtools view -bS -@ {} -'}
    ],
    ('bam', 'r'): [
        {'tool': 'samtools', 'command': 'samtools view -@ {} -h'}
    ],
    ('gz', 'w'): [
        {'tool': 'pbgzip', 'command': 'pbgzip -c -n {}'},
        {'tool': 'bgzip', 'command': 'bgzip -c -@ {}'},
        {'tool': 'gzip', 'command': 'gzip -c'}
    ],
    ('gz', 'a'): [
        {'tool': 'pbgzip', 'command': 'pbgzip -c -n {}'},
        {'tool': 'bgzip', 'command': 'bgzip -c -@ {}'},
        {'tool': 'gzip', 'command': 'gzip -c'}
    ],
    ('gz', 'r'): [
        {'tool': 'pbgzip', 'command': 'pbgzip -dc -n {}'},
        {'tool': 'bgzip', 'command': 'bgzip -dc -@ {}'},
        {'tool': 'gzip', 'command': 'gzip -dc'}
    ],
    ('lz4', 'w'): [
        {'tool': 'lz4c', 'command': 'lz4c -cz'}
    ],
    ('lz4', 'a'): [
        {'tool': 'lz4c', 'command': 'lz4c -cz'}
    ],
    ('lz4', 'r'): [
        {'tool': 'lz4c', 'command': 'lz4c -cd'}
    ]
}


@dataclass
class CommandRunResult:
    """
    CommandRunResult is a simple class that represents the IO command execution.
    It encapsulates the standard input, output, and error streams, along with the mode of operation.

    Mode of operation dictates the choice of input/output files: 
        'r' = reading, when output file is disabled,
        'w' = writing, when the input file is disabled,
        'a' = appending (input file is disabled).

    Attributes:
        input (Optional[TextIO]): The standard input stream (stdin) of the process. Can be None.
        errors (Optional[TextIO]): The standard error stream (stderr) of the process. Can be None.
        output (Optional[TextIO]): The standard output stream (stdout) of the process. Can be None.
        mode (Optional[Literal['r', 'w', 'a']]): The mode in which the file is opened: 'r' for reading, 'w' for writing, or 'a' for appending.

    Properties:
        outfile (Optional[TextIO]): Returns the appropriate file-like object based on the mode:
            - 'r': Returns the output stream (stdout) for reading.
            - 'w' or 'a': Returns the input stream (stdin) for writing or appending.
    """

    input: tp.Optional[tp.TextIO]
    errors: tp.Optional[tp.TextIO]
    output: tp.Optional[tp.TextIO]
    mode: tp.Optional[tp.Literal['r', 'w', 'a']]

    @property
    def outfile(self) -> tp.Optional[tp.TextIO]:
        if self.mode=='r':
            return self.output
        if self.mode=='w' or self.mode=='a':
            return self.input


@dataclass
class CommandFormatter():
    """
    CommandFormatter is a class that manages file opening operations with support for various compression formats.

    Attributes:
        mode (str): Mode in which the file is to be opened ('r', 'w', or 'a').
        path (Optional[str]): Path to the target file. Empty (None) or '-' indicates standard input/output (for the case of piping).
        command (Optional[Union[List[str], str]]): Custom command for file processing. For some file formats we have default commands. If empty or None, the class will try to find a suitable command based on the file format and mode.
        If a command is provided, it will be used directly.
        nproc (int): Number of threads for multithreaded tools. Defaults to 1.

    Methods:
        __call__(): Executes the command or opens the file based on the provided parameters. Return pairtools.lib.fileio.CommandRunResult object
    """

    mode: tp.Literal['r', 'w', 'a']
    path: tp.Optional[str]=None
    command: tp.Optional[tp.Union[tp.List[str], str]]=None
    nproc: int=1

    def __post_init__(self):
        """
        Auto-detect file extension, check and pick the tool available in the system.
        """
        self.__nocommand = False

        # If the command was provided by the user, simply run it, no further modifications needed:
        if self.command:
            return
        
        # If no user-defined command was provided, detect the command automatically.
        # Get the file extension:
        if not self.path or self.path == '-':
            self.__extension = ''
        else:
            self.__extension = self.path.split('.')[-1]

        # If extension is not in the keys of PRESET_EXT2MODE, 
        # simply return opened file in a given mode
        if self.__extension not in PRESET_EXT2MODE.keys():
            self.__nocommand = True
            return
        
        # If given mode is not applicable to the given extension, raise an error:
        if (self.__extension, self.mode) not in PRESET_COMMANDS.keys():
            raise ValueError(f'{self.__extension} can not to be opened in "{self.mode}" mode')
        
        # Next, iterate over possible tools. When the tool is found, constructs a command.
        checked_tools = []
        for possible_solution in PRESET_COMMANDS[(self.__extension, self.mode)]:
            # Check for the presence of the tool:
            if shutil.which(possible_solution['tool']) is None:
                checked_tools.append(possible_solution['tool'])
                continue

            # bgzip sometimes is outdated in the system, causing -@ option to crash. 
            # we explicitly check bgzip version, and if it is not up-to-date, skip bgzip:
            if possible_solution['tool']=='bgzip':
                cmd_check = subprocess.Popen(shlex.split('bgzip --version'), stderr=subprocess.PIPE, text=True)
                cmd_error = cmd_check.stderr.buffer.readline().decode()
                if 'invalid' in cmd_error:
                    checked_tools.append(possible_solution['tool'])
                    continue

            self.command = possible_solution['command'].format(str(self.nproc))
            return

        # No suitable command was found in the system, raise and format an error message.
        raise ValueError(self.format_notfounderror(checked_tools, self.mode=='r'))

    @staticmethod
    def format_notfounderror(checked_tools: tp.List[str], mode: bool) -> str:
        """
        Format a neat error message, used in case when none tools were found in the system.
        """
        text_task = 'read input file' if mode=='r' else 'write input file'
        if len(checked_tools)==0:
            raise ValueError('Something went wrong while IO operations. None tools were checked.')
        elif len(checked_tools)==1:
            text_verb = 'is'
            text_tools = checked_tools[0]
        else:
            text_verb = 'are'
            text_tools = f'{"", "".join(checked_tools[:-1])} and {checked_tools[-1]}'
        return f"{text_tools} {text_verb} not found, cannot {text_task}"
    
    def __construct_process(self):
        """
        Construct subprocess Popen object for a command, file path and file opening mode.
        """

        # Check if command is a string and contains shell operators like pipes
        use_shell = False
        if isinstance(self.command, str):
            if '|' in self.command or '>' in self.command or '<' in self.command:
                use_shell = True
                command_exec = self.command
            else:
                command_exec = shlex.split(self.command)
        else:
            # If already a list, assume no shell needed:
            command_exec = self.command

        # Open the file:
        self.__process_file = open(self.path, self.mode)

        # Run Popen:
        if self.mode == 'r':
            cmd = subprocess.Popen(command_exec, stdin=self.__process_file, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, shell=use_shell)
        else:
            cmd = subprocess.Popen(command_exec, stdout=self.__process_file, stdin=subprocess.PIPE, stderr=subprocess.PIPE, text=True, shell=use_shell)
        
        return cmd
    
    def __call__(self) -> CommandRunResult:
        """
        Execute the command.
        """

        # Empty path or path '-' means that we just read from stdin or write to stdout (in case of piping)
        if not self.path or self.path == "-":
            # sys.stdout and sys.stdin are inverted:
            return CommandRunResult(input=sys.stdout, errors=None, output=sys.stdin, mode=self.mode)

        # No command can happen when: 
        # (1) the provided extension is not in PRESET_EXT2MODE 
        # and (2) no command was provided by user.
        # We then simply open file in given mode and do not run any command on it:
        if self.__nocommand:
            self.__process_file = open(self.path, self.mode)
            return CommandRunResult(input=self.__process_file, errors=None, output=self.__process_file, mode=self.mode)

        process = self.__construct_process()
        return CommandRunResult(input=process.stdin, errors=process.stderr, output=process.stdout, mode=self.mode)


def auto_open(path, mode, nproc=1, command=None):
    """
    Automatically opens a file based on its format and access mode.

    Determines the file format from its extension and selects the appropriate
    command for reading or writing, utilizing available compression or
    decompression tools. If user defines an input command, it will be used
    instead.

    Parameters:
        path (str): Path to the file or '-' for standard input/output. Can be None in a piping regime.
        mode (str): File access mode: 'r' (read), 'w' (write), or 'a' (append).
        nproc (int, optional): Number of threads for multithreaded tools. Default is 1.
        command (str or list, optional): Custom command for file processing. Default is None.

    Returns:
        file-like object: A file object ready for reading or writing data.

    Raises:
        ValueError: File format detected from extension is not supported.
        ValueError: Required tool is not found available in the system.
    """
    command = CommandFormatter(mode=mode, path=path, command=command, nproc=nproc)
    result = command()
    return result.outfile


def get_stream_handlers(instream):
    """
    Get the readline and peek functions for the provided input stream.

    Parameters:
        instream (file-like object): The input stream to get the handlers for.

    Returns:
        tuple: A tuple containing the following elements:
            - readline_f (function): The readline function for the input stream.
            - peek_f (function): The peek function for the input stream.

    Raises:
        ValueError: If the peek function cannot be found for the provided stream.
    """
    readline_f, peek_f = None, None
    if hasattr(instream, "buffer"):
        peek_f = instream.buffer.peek
        readline_f = instream.buffer.readline
    elif hasattr(instream, "peek"):
        peek_f = instream.peek
        readline_f = instream.readline
    else:
        raise ValueError("Cannot find the peek() function of the provided stream!")
    return readline_f, peek_f


#### Legacy code (consider removing):

class PipedIO:
    def __init__(self, file_or_path, command, mode="r"):
        """
        An experimental class that reads/writes a file, piping the contents
        through another process.

        Parameters
        ----------
        file_or_path : file-like object or str
            A path to the input/output file or an already opened
            file-like object.
        command : str
            A command to launch a reading/writing process.
            If mode is 'w', the process must accept input via stdin.
            If mode is 'r', the process must put output into stdout.
            If mode is 'r' and file_or_path is str, the path will be
            appended to the command as the last argument.
        mode : str
            The mode for opening, same as in open(mode=).

        Returns
        -------
        file: a file-like object
        """

        if issubclass(type(command), str):
            command = command.split(" ")
        self._command = command
        self._mode = mode

        if mode.startswith("r"):
            if issubclass(type(file_or_path), str):
                self._proc = subprocess.Popen(
                    command + [file_or_path],
                    universal_newlines=True,
                    stdout=subprocess.PIPE,
                )
            else:
                self._proc = subprocess.Popen(
                    command,
                    universal_newlines=True,
                    stdin=file_or_path,
                    stdout=subprocess.PIPE,
                )
            self._stream = self._proc.stdout

            self._close_stream = self._proc.stdout.close

        elif mode.startswith("w") or mode.startswith("a"):
            f = (
                open(file_or_path, mode=mode)
                if issubclass(type(file_or_path), str)
                else file_or_path
            )
            self._proc = subprocess.Popen(
                command, universal_newlines=True, stdin=subprocess.PIPE, stdout=f
            )
            self._stream = self._proc.stdin

        self.buffer = self._stream.buffer
        self.closed = self._stream.closed
        self.flush = self._stream.flush
        self.fileno = self._stream.fileno

        self.read = self._stream.read
        self.readline = self._stream.readline
        self.readlines = self._stream.readlines

        self.seek = self._stream.seek
        self.seekable = self._stream.seekable
        self.truncate = self._stream.truncate
        self.tell = self._stream.tell

        self.writable = self._stream.writable
        self.write = self._stream.write
        self.writelines = self._stream.writelines

    def close(self, timeout=None):
        self._stream.close()
        retcode = self._proc.wait(timeout=timeout)
        return retcode
