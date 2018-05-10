#!/usr/bin/env python
import nibabel as nib
import numpy as np
import pandas as pd
from vtk.util.numpy_support import vtk_to_numpy
import vtk
from scipy.interpolate import interp1d
from tqdm import tqdm
import sys
import os

def importVTP(filename):
    streamlines = [] #list of streamlines
    filename = str(filename)
    vwriter = vtk.vtkXMLPolyDataWriter()
    reader = vtk.vtkXMLPolyDataReader()
    reader.SetFileName(filename)
    reader.Update()
    reader.ReleaseDataFlagOn()
    data = reader.GetOutput()
    vtxs = vtk_to_numpy(data.GetPoints().GetData())
    pts = data.GetPointData()
    vals = {pts.GetArrayName(idx):vtk_to_numpy(pts.GetArray(idx)) for idx in range(pts.GetNumberOfArrays())} # dictionary of numpy array of scalar values associated with each vertex
    metrics = [pts.GetArrayName(idx) for idx in range(pts.GetNumberOfArrays())] # name of metrics
    if 'tensors' in metrics:
        del metrics[metrics.index('tensors')]
    for i in range(data.GetNumberOfCells()):
        streamlines.append([data.GetCell(i).GetPointIds().GetId(p) for p in range(data.GetCell(i).GetPointIds().GetNumberOfIds())])
    return(vtxs,vals,streamlines)

def saveVTP(vtxs,vals,streamlines,filename):
    polydata = vtk.vtkPolyData()
    points = vtk.vtkPoints()
    lines = vtk.vtkCellArray()
    points.SetNumberOfPoints(len(vtxs))
    vtxs2 = vtxs
    for i, p in tqdm(enumerate(vtxs2)):
        points.SetPoint(i,p[0],p[1],p[2])
    for stream in streamlines:
        lines.InsertNextCell(len(stream))
        for i in stream:
            lines.InsertCellPoint(i)
    polydata.SetPoints(points)
    polydata.SetLines(lines)
    pointdata = polydata.GetPointData()
    for sname, sarr in vals.iteritems():
        arr = vtk.vtkFloatArray()
        arr.SetName(sname)
        arr.SetNumberOfComponents(1)
        for v in sarr:
            arr.InsertNextTuple1(v)
        pointdata.AddArray(arr)
        pointdata.SetActiveScalars(sname)
    vwriter = vtk.vtkXMLPolyDataWriter()
    vwriter.SetInput(polydata)
    vwriter.SetFileName(str(filename))
    vwriter.Write()

def point_affine(arr,T1file):
    loli = nb.Nifti1Image.load(T1file)
    loli_affine = loli.affine
    pts = nb.affines.apply_affine(np.linalg.inv(loli_affine),np.array(arr))
    return(pts)

def extract_from_ROI(ROI_file,vtxs,vals,**diff_dic):
    global use_MM
    ROI_file = str(ROI_file)
    data = nib.load(ROI_file)
    img = data.get_data()

    # collate points within ROI
    if use_MM:
        ROI_points = []
        ROI_values = {'AD':[],'RD':[],'MD':[],'FA':[], 'MM':[]}
        rest_points = []
        rest_values = {'AD':[],'RD':[],'MD':[],'FA':[], 'MM':[]}
    else:
        ROI_points = []
        ROI_values = {'AD':[],'RD':[],'MD':[],'FA':[]}
        rest_points = []
        rest_values = {'AD':[],'RD':[],'MD':[],'FA':[]}

    voxs = np.vstack(np.where(img == 1)).T #get all voxels labeled 1
    count = 0
    for vox in voxs:
        maxs = np.apply_along_axis(lambda x: x+0.5,0,vox)
        mins = np.apply_along_axis(lambda x: x-0.5,0,vox)
        maxs = nib.affines.apply_affine(data.affine,maxs)
        mins = nib.affines.apply_affine(data.affine,mins)
        mins,maxs = MinMax_reorganization(mins,maxs)
        count+=1
        for i in tqdm(range(len(vtxs)),desc='Processing voxel '+str(count)):
            if mins[0] <= vtxs[i][0] < maxs[0]:
                if mins[1] <= vtxs[i][1] < maxs[1]:
                    if mins[2] <= vtxs[i][2] < maxs[2]:
                        ROI_points.append(vtxs[i])
                        for k in ROI_values.keys():
                            ROI_values[k].append(vals[k][i])
                            #print(vals[k][i])
            else:
                rest_points.append(vtxs[i])
                for k in rest_values.keys():
                    rest_values[k].append(vals[k][i])
    ROI_points = np.array(ROI_points)
    rest_points = np.array(rest_points)
    for k in ROI_values.keys():
        ROI_values[k] = np.array(ROI_values[k])
    for k in rest_values.keys():
        rest_values[k] = np.array(rest_values[k])

    return(ROI_points,ROI_values,rest_points,rest_values)


def extract_from_sphere(ROI_file,vtxs,vals,radius=1.5):
    global use_MM
    ROI_file = str(ROI_file)
    data = nib.load(ROI_file)
    img = data.get_data()
    means = np.mean(np.vstack(np.where(img == 1)),axis=1) # obtain the center of the ROI_file
    #apply affines
    means = nib.affines.apply_affine(data.affine,means)
    # collate points within sphere
    if use_MM:
        ROI_points = []
        ROI_values = {'AD':[],'RD':[],'MD':[],'FA':[], 'MM':[]}
        rest_points = []
        rest_values = {'AD':[],'RD':[],'MD':[],'FA':[], 'MM':[]}
    else:
        ROI_points = []
        ROI_values = {'AD':[],'RD':[],'MD':[],'FA':[]}
        rest_points = []
        rest_values = {'AD':[],'RD':[],'MD':[],'FA':[]}
    for i in tqdm(range(len(vtxs))):
        if np.abs(np.linalg.norm(vtxs[i]-means)) <= radius:
            ROI_points.append(vtxs[i])
            for k in ROI_values.keys():
                ROI_values[k].append(vals[k][i])
        else:
            rest_points.append(vtxs[i])
            for k in rest_values.keys():
                rest_values[k].append(vals[k][i])
    ROI_points = np.array(ROI_points)
    rest_points = np.array(rest_points)
    for k in ROI_values.keys():
        ROI_values[k] = np.array(ROI_values[k])
    for k in rest_values.keys():
        rest_values[k] = np.array(rest_values[k])
    return(ROI_points,ROI_values,rest_points,rest_values)

def resample_data(vals,timep,n=1000,diffs=['AD','RD','MD','FA']):
    df=pd.DataFrame()
    timep = str(timep)
    print('# of points: %s' % (vals[diffs[0]].shape[0]))
    if vals[diffs[0]].shape[0] < n:
        n = vals[diffs[0]].shape[0]
        print('Downsampling capped at: ' +str(n))
    elif n == False:
        n = vals[diffs[0]].shape[0]
        print('Not downsampling.')
    for diff in diffs:
        if n:
            y = np.sort(vals[diff])
            x = np.arange(0,y.shape[0],1)
            f = interp1d(x,y)
            values = f(np.linspace(0,y.shape[0],n,endpoint=False))
            #values = resample(np.sort(vals[diff]),n)
        else:
            values = vals[diff]
        df[diff] = values
        df['Time'] = np.array([timep for e in range(len(values))])
    df = df.drop('Time',1)
    return(df)

def MinMax_reorganization(mins,maxs):
    if maxs[0] > mins[0]:
        XMax = maxs[0]
        XMin = mins[0]
    else:
        XMax = mins[0]
        XMin = maxs[0]
    if maxs[1] > mins[1]:
        YMax = maxs[1]
        YMin = mins[1]
    else:
        YMax = mins[1]
        YMin = maxs[1]
    if maxs[2] > mins[2]:
        ZMax = maxs[2]
        ZMin = mins[2]
    else:
        ZMax = mins[2]
        ZMin = maxs[2]
    return(np.array([XMin,YMin,ZMin]),np.array([XMax,YMax,ZMax]))
'''Start main process'''

#pymp.config.nested=True
try:
    input_tractname = str(sys.argv[1])
    input_roi1 = str(sys.argv[2])
    output_filename = str(sys.argv[3])
    output_prefix = output_filename.replace('.csv','')+'-'
    use_sphere = int(sys.argv[4])
    if use_sphere == 0:
        use_sphere = False
    use_MM = int(sys.argv[5])
    if use_MM == 0:
        use_MM = False
    vtxs,vals,_=importVTP(input_tractname)
    print(vals['AD'].shape)
    vtxs2=vtxs
    vals2=valscut_plane.actor.property.opacity = 0.7
except:
    print('vtxStats.py | Usage: <input_tractname> <input_roi> <output_filename> <use_sphere> <use_MM> [n_samples=1000] [radius=2mm]')
    print('Version May 10 2018\nextract_from_ROI modified to extract from every labeled voxels independently.')
    #print('Version Feb 20 2018\nUse_MM option added.')
    #print('Version Mar 15 2018\nCSV output option')
    sys.exit(2)
if use_MM:
    diff_list = ['AD','RD','MD','FA','MM']
else:
    diff_list = ['AD','RD','MD','FA']
try:
    n_samples = int(sys.argv[6])
    if n_samples == 0:
        n_samples = False
except:
    n_samples = 1000
try:
    radius = float(sys.argv[7])
except:
    radius = 2
'''print(sys.argv[1])
print(sys.argv[2])
print(sys.argv[3])
print(sys.argv[4])
print(sys.argv[5])
print(sys.argv[6])
print(sys.argv[7])'''

os.system('mkdir -p stats')

if use_sphere:
    _,values,_,_=extract_from_sphere(input_roi1,vtxs,vals,radius=radius)
    print('Spherical extraction')
else:
    _,values,_,_=extract_from_ROI(input_roi1,vtxs,vals)
df = resample_data(values,timep="",n=n_samples,diffs=diff_list)
for diff in diff_list:
    tract_mean = np.nanmean(df[diff])
    print("%s: %s" % (diff,df[diff].mean()))
    if tract_mean != np.nan:
        cur_subj = input_tractname.replace('.vtp','')
        os.system('echo '+cur_subj+', '+str(tract_mean)+' >> stats/'+output_prefix+diff+'.csv')
