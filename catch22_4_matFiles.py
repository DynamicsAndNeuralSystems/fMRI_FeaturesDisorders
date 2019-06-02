import numpy as np
import python2Functions as p2f

# Assign a path (make sure you add that last slash) and name the dataset
path = '/Users/AV/Dropbox/COBRE/cfgData/'
dataSet = 'COBRE'

# Select which element of 'roiTS' should be analysed by changing the second input
x, element = p2f.giveMeBigFeatMat(path,3)

# Save file
np.savetxt('element' + str(element) + '_' + str(dataSet) + '.txt',x,delimiter=',')
