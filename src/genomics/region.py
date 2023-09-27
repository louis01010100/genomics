class Region():

    def __init__(self, start, end):
        self.start = start
        self.end = end

    def is_overlapping(self, other: Region) -> bool:
        return self.start <= other.end and other.start <= self.end
