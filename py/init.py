import sys
import pyR2D2

# Initialize R2D2.Data object
if __name__ == '__main__':
    if 'google' in sys.argv:
        google = True
    else:
        google = False
    
    pyR2D2.util.init(locals(),google=google)