import time


class Timer(object):
    def __init__(self, name=None):
        self.name = name

    def __enter__(self):
        self.tstart = time.time()

    def __exit__(self, type, value, traceback):
        True
    #     if self.name:
    #         print('[%s]' % self.name,)
    #     print('Elapsed: %.2f seconds' % (time.time() - self.tstart))