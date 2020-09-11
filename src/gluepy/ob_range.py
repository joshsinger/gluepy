
class OneBasedRangeException(Exception):
    pass

class OneBasedRange:
    def __init__(self, start, end):
        if(start < 1):
            raise OneBasedRangeException("OneBasedRange cannot have negative start")
        if(end < start):
            raise OneBasedRangeException("OneBasedRange start may not be greater than end")
        self.start = start
        self.end = end
    def __str__(self):
        return "ob["+str(self.start)+":"+str(self.end)+"]"
    def pull_from_str(self, string):
        if(self.end > len(string)):
            raise OneBasedRangeException("Cannot use OneBasedRange "+str(self)+" on FastaSequence of length "+str(len(string)))
        return string[(self.start-1):(self.end)]
    def length(self):
        return (self.end - self.start)+1