"""Programmatic access to ClarioStar plate reader data"""
import re
import numpy
import xlrd


class Assay:
    """Plate reader assay

    Takes an excel file from ClarioStar (XXX Model?) and gives read access
    to the contained data:

    >>> assay = Assay(path_to_file, split=True)

    The Assay instance will have the following attributes:

    assay.time,  a numpy array of recorded times
    assay.well,  an array with fliorescence measurements for each well
    assay.phase, list of time ranges with contiguous recordings

    Wells are adressed by a single integer index ranging from 0 to 95, e.g.

    >>> assay.well[14]

    gives the timeline for well 14 (colum B row 2 in a 12x8 layout).

    If split is True (the default), then the raw data is grouped into
    time ranges of contiguous recordings. assay.phase[0] contains the entire
    time interval, and assay.phase[1:] contains ranges of individual contigs.
    """
    COLS = 12
    ROWS = 8

    def __init__(self, path, split=True):
        def parse_time(string):
            match = re.match(r'Cycle (.*) \((.*) min( (.*) s)?\)', string)
            mins = int(match.group(2))
            secs = int(match.group(4)) if match.group(4) else 0
            return mins + secs/60

        book = xlrd.open_workbook(path)
        ptcl = book.sheet_by_index(2)
        data = book.sheet_by_index(0)

        # read protocol
        self.test_name = ptcl.cell(3, 0).value[len('Test Name: '):]
        self.measurement = ptcl.cell(11, 1).value
        self.cycles = int(ptcl.cell(17, 1).value)
        self.cycle_time = float(ptcl.cell(18, 1).value/60)

        # initialize data attributes
        self.time = numpy.zeros(self.cycles)
        self.well = numpy.zeros((self.ROWS*self.COLS, self.cycles))

        for idx in range(self.cycles):
            # parse time
            self.time[idx] = parse_time(data.cell(12 + (self.ROWS+4)*idx, 0).value)

            # read wells
            for col in range(self.COLS):
                for row in range(self.ROWS):
                    val = data.cell(15 + (self.ROWS+4)*idx + row,
                                    1 + col).value
                    self.well[row*self.COLS + col, idx] = val if val != '' else None

        # split time into consecutive phases
        self.phase = [range(0, len(self.time))]
        if split:
            # XXX deletes the last data point of every phase
            mask = numpy.diff(self.time) > self.cycle_time + 0.1/60
            mask = numpy.append(mask, False)
            self.time = numpy.ma.array(self.time)
            self.time[mask] = numpy.ma.masked

            left = 0
            size = -1
            while size:
                size = self.time.mask[left:].argmax()
                self.phase.append(range(left, left+size if size else len(self.time)))
                left += size + 1
