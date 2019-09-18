import numpy as np
import python2Functions as p2f

# Assign a path (make sure you add that last slash) and name the dataset
path = '/Users/AV/Dropbox/UCLA/cfgData/'
dataset = 'UCLA'

# Select which element of 'roiTS' should be analysed by changing the second input
x, element = p2f.giveMeBigFeatMat(path,1)

# Save file
np.savetxt('element' + str(element) + '_' + str(dataset) + '.txt',x,delimiter=',')
