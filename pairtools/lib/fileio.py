import shutil
import subprocess
import shlex
import sys
import typing as tp
from dataclasses import dataclass

class ParseError(Exception):
    pass


MODES_TO_FILES_PRESET = {
    'bam': ['w', 'r'],
    'gz': ['w', 'r', 'a'],
    'lz4': ['w', 'r', 'a'],
}

COMMANDS = {
    ('bam', 'w'): [
        {'tool': 'samtools', 'command': 'samtools view -bS {} -'}
    ],
    ('bam', 'r'): [
        {'tool': 'samtools', 'command': 'samtools view -h'}
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
    errors: tp.Optional[bytes]
    output: tp.Optional[bytes]
    outfile: tp.Optional[tp.TextIO]


@dataclass
class CommandFormatter():

    mode: tp.Optional[tp.Literal['r', 'w', 'a']]
    path: tp.Optional[str]=None
    command: tp.Optional[tp.Union[tp.List[str], str]]=None
    nproc: int=1
    is_binary: bool=False

    @staticmethod
    def form_notfounderror_text(searched_tools: tp.List[str], is_read: bool) -> str:
        tools_article = 'is' if len(searched_tools) == 1 else 'are'
        tools_defenition = 'compress output' if is_read else 'decompress input'
        tools_list = f'{"", "".join(searched_tools[:-1])} and {searched_tools[-1]}' if len(searched_tools) > 1 else searched_tools[0]
        return f"{tools_list} {tools_article} not found, cannot {tools_defenition}"

    def __post_init__(self):
        self.__nocommand = False
        if self.is_binary:
            self.file_mode = f'{self.mode}b'
        else:
            self.file_mode = self.mode

        if self.command:
            return
        
        self.__format = self.path.split('.')[-1]

        if self.__format not in MODES_TO_FILES_PRESET.keys():
            self.__nocommand = True
            return
        if (self.__format, self.mode) not in COMMANDS.keys():
            raise ValueError(f'{self.__format} can not to be opened in {self.mode}')
        
        checked_tools = []
        for possible_solution in COMMANDS[(self.__format, self.mode)]:
            if shutil.which(possible_solution['tool']) is None:
                checked_tools.append(possible_solution['tool'])
                continue
            self.command = possible_solution['command'].format(str(self.nproc))
            return

        raise ValueError(self.form_notfounderror_text(checked_tools, self.mode=='r'))
    
    def __convert_command(self):
        if isinstance(self.command, str):
            self.__command_to_sp = shlex.split(self.command)
        else:
            self.__command_to_sp = self.command

    def __form_command(self):
        self.__process_file = open(self.path, self.file_mode)
        if self.mode == 'r':
            print('read')
            cmd=subprocess.Popen(self.__command_to_sp, stdin=self.__process_file, stdout=subprocess.PIPE)
        else:
            print('write')
            cmd=subprocess.Popen(self.__command_to_sp, stdout=self.__process_file, stdin=subprocess.PIPE)
        return cmd
    
    def __call__(self):
        if not self.path or self.path == "-":
            if self.mode == "r":
                return CommandRunResult(errors=None, output=None, outfile=sys.stdin)
            if self.mode == "w":
                return CommandRunResult(errors=None, output=None, outfile=sys.stdout)
            if self.mode == "r":
                return ValueError(f'Zero input file can not to be opened in a')
        if self.__nocommand:
            return CommandRunResult(errors=None, output=None, outfile=open(self.path, self.file_mode))
        self.__convert_command()
        process = self.__form_command()
        out, err = process.communicate()
        return CommandRunResult(errors=err, output=out, outfile=self.__process_file)


def auto_open(path, mode, nproc=1, command=None):
    command = CommandFormatter(mode=mode, path=path, command=command, nproc=nproc)
    result = command()
    return result.outfile


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


