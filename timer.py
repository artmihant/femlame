from datetime import datetime
class Timer:
    def __init__(self,message):
        self.now = now()
        self.start = now()
        print('\t',message)
    def __call__(self,message):
        current = now()
        print(
            # time_print(self.now,current),'\t',
            time_print(self.start,current),'\t',message )
        self.now = current

def time_print(t1,t2):
    return '{:.1f}'.format((t2-t1).total_seconds())

def now():
    return datetime.now()
