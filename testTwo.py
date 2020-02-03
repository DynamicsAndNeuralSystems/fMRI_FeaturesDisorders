import pandas as pd
import numpy as np
from scipy import stats
from statsmodels.stats import multitest
import sys
def pandasIndexTest():
    df = pd.DataFrame({'Yo':[1,1,1,3],'Nah':['a', 'b', 'c','d'], 'Yeah':[8,8, 9, 7]})
    df.set_index(['Yo'],inplace=True)
    nums = np.asarray(df.loc[1]['Yeah'])
    print(np.mean(nums>8))

def accuracyPValueTest():
    randomLearnDataUCLA = pd.read_csv("randomLearnData_UCLA.txt")
    randomLearnDataCOBRE = pd.read_csv("randomLearnData_COBRE.txt")
    # featAccs = np.asarray([1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22])
    randomLearnDataUCLA.set_index('Feature',inplace=True)
    randomLearnDataCOBRE.set_index('Feature',inplace=True)
    mockPValUCLA = np.zeros(22)
    mockPValCOBRE = np.zeros(22)
    for feature in range(22):
        randomLearnSliceUCLA = randomLearnDataUCLA.loc[feature]
        randomLearnSliceCOBRE = randomLearnDataCOBRE.loc[feature]
        randomAccsUCLA = np.asarray(randomLearnSliceUCLA['Average Accuracy'])
        randomAccsCOBRE = np.asarray(randomLearnSliceCOBRE['Average Accuracy'])
        mockPValUCLA[feature] = np.mean(randomAccsUCLA>=50)
        mockPValCOBRE[feature] = np.mean(randomAccsCOBRE>=50)

    mockPValUCLA = multitest.multipletests(mockPValUCLA,method='fdr_bh')[1]
    mockPValCOBRE = multitest.multipletests(mockPValCOBRE,method='fdr_bh')[1]
    pVals = np.column_stack((mockPValUCLA,mockPValCOBRE))
    print(pVals)
    for pVal in pVals:
        smth = "%.3f" % stats.combine_pvalues(pVal, method='fisher')[1]
        print(smth)

def readProcMeth3LabelColumnTest():
    labels = np.genfromtxt('labels_procMeth3_UCLA.txt')
    c22NDArray = np.genfromtxt('c22_procMeth3_UCLA.txt')
    # fdAvgs = '/Users/preethompal/Dropbox (Sydney Uni Student)/UCLA/movementData/fdAvgs_UCLA.txt'
    # fdAvgs = pd.read_csv(fdAvgs, names=['Averages'])
    # print(fdAvgs)
    print((len(c22NDArray)/82))
    print('ctrls', np.sum(labels == 0))
    print('sczs', np.sum(labels == 1))
    ls = [i for i in range(int(len(c22NDArray)/82))]
    print(ls)

def numpyColumnTest():
    a = np.asarray([[True,True],[True,False],[True,False]])
    print(a[:,1].mean())

def dropMultiIndexRow():
    iterables = [['bar', 'baz', 'foo', 'qux'], list(range(10))]
    mdx = pd.MultiIndex.from_product(iterables, names=['first', 'second'])
    d = np.zeros(400).reshape((40,10))
    df = pd.DataFrame(data=d, index=mdx)
    print(df)
    indsToDrop = [1,2,3,4,5]
    df.drop(index=indsToDrop,level='second',inplace=True)

    index = df.index
    indices = []
    indices.append(list(index.levels[0]))
    indices.append([i for i in range(5)])
    idx = pd.MultiIndex.from_product(indices, names=['first', 'second'])
    df.index = idx
    # tdf = pd.DataFrame(data=df[i for i in df.columns], index=idx)
    print(df)

def dropNumpyArrayRow():
    a = np.zeros(10)
    a = a.reshape((len(a),1))
    print(len(a))
    print(a)
    niggers = [0,1,2,3]
    a = np.delete(a,niggers)
    a = a.reshape((len(a),1))
    print(len(a))
    print(a)

def dataframeColumnIterationTest():
    df = pd.DataFrame({'a':[float('NaN'), 1],'b':[2,2],'c':[3,2]})
    for featName, featColumn in df.iteritems():
        if featColumn.isnull().values.any():
            print('Null in column found.')
            print("Presence of NaNs: ", str(np.mean(featColumn.isnull().values)*100)+'%')
            df.drop(featName, axis='columns', inplace=True)
        else:
            print("Null-free column found.")
    print(df)

def indexOnFeatSliceShapeTest():

    df = pd.DataFrame({'Reg 1':[1,2,3], 'Reg 2':[4,5,6], 3:[7,8,9]})
    train_indices = [1,2]
    test_indices = [3]
    print(df[3])

def catch22AttributesTest():
    import catch22
    features = dir(catch22)
    features = [feat for feat in features if not '__' in feat]
    del features[22:24]
    print(features)

def numpyVstackTest():
    a = np.asarray([1,2,3,4])
    b = np.asarray([5,6,7,8])
    c = np.vstack([a,b])
    d = np.column_stack([a,b])
    print(c)
    print(d)

def multiIndexDFSliceTest():
    iterables = [['bar', 'baz', 'foo', 'qux'], list(range(10))]
    mdx = pd.MultiIndex.from_product(iterables, names=['regions', 'subjects'])
    d = np.zeros(400).reshape((40,10))
    newCol = list(range(40))
    df = pd.DataFrame(data=d, index=mdx)
    df[2] = newCol
    print(df)

    slice = df.loc[:,2].values.reshape(len(iterables[0]), len(iterables[1])).transpose()


    print(slice)
    df = pd.DataFrame(data=slice, columns = iterables[0])
    df.index.name='Subjects'
    a = np.asarray([[0],[0],[0],[0], [1],[1],[1],[1],[1],[1]])
    ctrls = np.where(a==0)[0]
    print(ctrls)
    print(df.iloc[ctrls, [1]])

def labelColumnIndicesSplitTest():
    a = np.asarray([[0],[0],[0],[0],[0],[1],[1],[1],[1],[1],])
    ctrls = np.where(a==0)[0]
    print(ctrls)

def toCSVAppendTest():
    df = pd.DataFrame({'Reg 1':[1,2,3], 'Reg 2':[4,5,6], 'Reg 3':[7,8,9]})
    outFileName='yo.txt'
    df.to_csv(outFileName, index=False)

    df1 = pd.DataFrame({'Reg 1':[1,2,3], 'Reg 2':[4,5,6], 'Reg 3':[7,8,9]})
    df1.to_csv(outFileName, mode='a',header=False,index=False)

def makeSuperMatrixTest():
    np.set_printoptions(threshold=sys.maxsize)
    a = np.zeros((5, 10*22))
    print(a.shape)

    iterables = [list(range(10)), list(range(5))]
    mdx = pd.MultiIndex.from_product(iterables, names=['regions', 'subjects'])
    d = np.zeros((10*5,22))
    df = pd.DataFrame(data=d, index=mdx)

    # for j in range(5):
    #     for i in range(10):
    #         a[i,22*j:22*(j+1)] = j

    for j in range(5):
        k = 0
        for i in range(10):
            df.loc[i,j] = j + k
            k+=1
            a[j,i*22:22*(i+1)] = df.loc[i,j]

    print(df)
    df = pd.DataFrame(data=a)
    df.index.name='Subjects'
    df.rename_axis('Regions,Features', axis='columns',inplace=True)
    print(df)

def numpyShuffleTest():
    a = np.asarray(list(range(10)))
    np.random.shuffle(a)
    print(a)
    ls = list(a)
    print(a.mean())


numpyShuffleTest()
