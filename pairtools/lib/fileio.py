import shutil
import pipes
import subprocess
import sys

class ParseError(Exception):
    pass


def auto_open(path, mode, nproc=1, command=None):
    """Guess the file format from the extension and use the corresponding binary
    to open it for reading or writing. If the extension is not known, open the
    file as text.

    If the binary allows parallel execution, specify the number of threads
    with `nproc`.

    If `command` is supplied, use it to open the file instead of auto-guessing.
    The command must accept the filename as the last argument, accept input
    through stdin and print output into stdout.

    Supported extensions and binaries (with comments):
    .bam - samtools view (allows parallel writing)
    .gz - pbgzip if available, otherwise bgzip
    .lz4 - lz4c (does not support parallel execution)
    """

    # Empty filepath or False provided
    if not path or path == "-":
        if mode == "r":
            return sys.stdin
        if mode == "w":
            return sys.stdout

    if command:
        if mode == "w":
            t = pipes.Template()
            t.append(command, "--")
            f = t.open(path, "w")
        elif mode == "r":
            t = pipes.Template()
            t.append(command, "--")
            f = t.open(path, "r")
        else:
            raise ValueError("Unknown mode : {}".format(mode))
        return f

    elif path.endswith(".bam"):
        if shutil.which("samtools") is None:
            raise ValueError(
                {
                    "w": "samtools is not found, cannot compress output",
                    "r": "samtools is not found, cannot decompress input",
                }[mode]
            )
        if mode == "w":
            t = pipes.Template()
            t.append(
                "samtools view -bS {} -".format(
                    "-@ " + str(nproc - 1) if nproc > 1 else ""
                ),
                "--",
            )
            f = t.open(path, "w")
        elif mode == "r":
            t = pipes.Template()
            t.append("samtools view -h", "--")
            f = t.open(path, "r")
        else:
            raise ValueError("Unknown mode for .bam : {}".format(mode))
        return f

    elif path.endswith(".gz"):
        if shutil.which("pbgzip") is not None:
            if mode == "w":
                t = pipes.Template()
                t.append("pbgzip -c -n {}".format(nproc), "--")
                f = t.open(path, "w")
            elif mode == "a":
                t = pipes.Template()
                t.append("pbgzip -c -n {} $IN >> $OUT".format(nproc), "ff")
                f = t.open(path, "w")
            elif mode == "r":
                t = pipes.Template()
                t.append("pbgzip -dc -n {}".format(nproc), "--")
                f = t.open(path, "r")
            else:
                raise ValueError("Unknown mode for .gz : {}".format(mode))
        elif shutil.which("bgzip") is not None:
            if mode == "w":
                t = pipes.Template()
                t.append("bgzip -c -@ {}".format(nproc), "--")
                f = t.open(path, "w")
            elif mode == "a":
                t = pipes.Template()
                t.append("bgzip -c -@ {} $IN >> $OUT".format(nproc), "ff")
                f = t.open(path, "w")
            elif mode == "r":
                t = pipes.Template()
                t.append("bgzip -dc -@ {}".format(nproc), "--")
                f = t.open(path, "r")
            else:
                raise ValueError("Unknown mode for .gz : {}".format(mode))
        elif shutil.which("gzip") is not None:
            if mode == "w":
                t = pipes.Template()
                t.append("gzip -c", "--")
                f = t.open(path, "w")
            elif mode == "a":
                t = pipes.Template()
                t.append("gzip -c $IN >> $OUT", "ff")
                f = t.open(path, "w")
            elif mode == "r":
                t = pipes.Template()
                t.append("gzip -dc", "--")
                f = t.open(path, "r")
            else:
                raise ValueError("Unknown mode for .gz : {}".format(mode))
        else:
            raise ValueError(
                {
                    "w": "pbgzip, bgzip and gzip are not found, cannot compress output",
                    "a": "pbgzip, bgzip and gzip are is not found, cannot compress output",
                    "r": "pbgzip, bgzip and gzip are is not found, cannot decompress input",
                }[mode]
            )
        return f
    elif path.endswith(".lz4"):
        if shutil.which("lz4c") is None:
            raise ValueError(
                {
                    "w": "lz4c is not found, cannot compress output",
                    "a": "lz4c is not found, cannot compress output",
                    "r": "lz4c is not found, cannot decompress input",
                }[mode]
            )
        if mode == "w":
            t = pipes.Template()
            t.append("lz4c -cz", "--")
            f = t.open(path, "w")
        elif mode == "a":
            t = pipes.Template()
            t.append("lz4c -cz $IN >> $OUT", "ff")
            f = t.open(path, "w")
        elif mode == "r":
            t = pipes.Template()
            t.append("lz4c -cd", "--")
            f = t.open(path, "r")
        else:
            raise ValueError("Unknown mode : {}".format(mode))
        return f
    else:
        return open(path, mode)


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


