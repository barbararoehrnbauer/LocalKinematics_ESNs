# -*- coding: mbcs -*-
from part import *
from material import *
from section import *
from assembly import *
from step import *
from interaction import *
from load import *
from mesh import *
from optimization import *
from job import *
from sketch import *
from visualization import *
from abaqus import *
from abaqusConstants import *
import customKernel #stores variables within the cernal itself 
import math
import random
import datetime
import operator



n_odb = 1 # Number of ESN to be generated per ESN configuration
b_ns = [40] # List of widths to be generated
t_ns = [1] # List of thickness to be generated
for t_i in t_ns:
    for b_i in b_ns:
        for i in range(n_odb):
            Mdb() # create new model to avoid errors if the script is executed multiple times
            ####################################################################################################################
            # parameter settings
            ###############################################################
            # model size parameters
            b_n = b_i# size of membrane patch (width & height) in um
            t_n = t_i  # thickness of membrane
            # crosslinking
            l_s_mean = 5 # mean segment length between two crosslinks
            # l_s_mean = False # Do not comment this out if you do not want to use l_mean in your calculations
            # loading parameters
            eps11 = 0.04
            eps12 = 0.04
            # choose one of three:
            LoadCase = 'Uniaxial'
            # LoadCase = 'PlanarTension'
            # LoadCase = 'SimpleShear'

            # choose one of two:   
            # #sinusoidal fibre generation with mirroring operations for PBC
            #FibreType = 'Sinusoidal' 
            # randomwalk fibre generation where PBC is linked between end of initial fibre & start of following fibre
            FibreType = 'RandomWalk'    

            ###############################################################
            # set parameters by known material
            # general
            # Material settings
            youngs = 1575 
            poisson = 0.3 
            plasticity = ((31.50, 0.0), (46.29, 0.017), (54.36, 0.0322), (59.44, 0.0497), (63.28, 0.0738), (66.24, 0.1198), (70.88,0.2031), (144,1.56)) #https://classes.engineering.wustl.edu/2009/spring/mase5513/abaqus/docs/v6.5/books/gsx/default.htm?startat=ch05s02.html
            d_f = 0.403 # Average fibre diameter measured from SEM images
            por = 0.87 # porosity of material,
            # curvature stuff
            # Isotropy
            Alg = 1     # Alignment value (0,1]: 1= isotropic, 0= all fibres left to right in one line
                        # only work in RW so far
            # sinusoidal specific
            xi_f_av = 7.42          # average amplitude of sine wave function describing a fibre in um
            xi_f_stdev = 4.82       # Standard deviation of amplitude
            delta_f_av = 130.57     # wavelength of sine wave function describing a fibre in um
            delta_f_stdev = 109.49  # Standard deviation of wavelength
            delta_f_min = 32        # minimum wavelength
            # RandomWalk specific
            ang_av = 0.138380645
            ang_stdev = 0.124080001
            ang_min = 0.001207973
            ang_max = 0.354377867

            ###############################################################
            # Performance parameters
            SplinePoints = 16 # nunmer of interpolation points within one wavelength of a fibre (sinusoidal)
            l_step = 3              # segment length in um (distance between two points of the spline)
            res_fp = 0.5 # fibre crosslinking range & subsequent partitioning size

            ###############################################################
            # calculated parameters
            TotalFibreLength = (1-por)*4*b_n**2*t_n/(math.pi*d_f**2) #choose this if you want to create a reality inspired model
            # TotalFibreLength = 2500 #choose this for debugging
            #Int_t = int(min([Int_t,TotalFibreLength/b_n*2-1])) #
            b_2 = float(b_n)/2
            # today's date
            date = datetime.date.today().strftime('%Y-%m-%d')
            counter = 1

            # save stuff in model file for later retrieval if necessary
            mdb.customData.savefloats = [b_n, eps11, eps12, res_fp, t_n,l_s_mean]
            mdb.customData.saveStrings = [FibreType,LoadCase]


            ######################################################################################################################



            # single creation steps:
            # model creation
            m = mdb.models['Model-1']

            # Material creation
            m.Material(name='fibreMat')
            m.materials['fibreMat'].Elastic(table=((youngs, poisson), ))
            m.materials['fibreMat'].Plastic(scaleStress=None, table=plasticity)
            m.materials['fibreMat'].Density(table=((1.78e-09, ), ))
            # section creation
            m.CircularProfile(name='Circular', r=d_f/2)
            m.BeamSection(consistentMassMatrix=False, integration=
                DURING_ANALYSIS, material='fibreMat', name='Circular_Beam', 
                poissonRatio=poisson, profile='Circular', temperatureVar=LINEAR)

            # step creation
            m.StaticStep(adaptiveDampingRatio=None, continueDampingFactors=False, 
                initialInc=0.02, maxInc=0.05, maxNumInc=10000, name='LoadingStep-1', 
                nlgeom=ON, previous='Initial', stabilizationMagnitude=0.0005, 
                stabilizationMethod=DISSIPATED_ENERGY_FRACTION)
            m.steps['LoadingStep-1'].setValues(minInc=1e-08)

            m.fieldOutputRequests['F-Output-1'].setValues(timeInterval=0.05,variables=(
                'S', 'MISES','PE', 'PEEQ', 'PEMAG', 'LE', 'U', 'RF', 'RT', 'CF', 'VF', 'BF', 
                'CSTRESS', 'CDISP', 'CFORCE', 'ENER', 'ELEN', 'ELEDEN','EE','LE','E','COORD'))

            m.historyOutputRequests['H-Output-1'].setValues(
                timeInterval=0.05, variables=('ALLSE', 'ALLSD', 'ETOTAL','ALLKE','ALLWK','ALLPD'))

            # Sinusoidal fibre generation
            if FibreType == 'Sinusoidal':
                # fibre creation
                AccLength = 0
                fibrecount = 0
                x_rand = []
                y_rand = []
                angle_rand = []
                xi_f =[]
                delta_f = []
                res_f = []
                while AccLength < TotalFibreLength: 
                    NoP = 'F'+str(fibrecount*4) # construct name of part

                    # tbd randomized parameters
                    x_rand.append(random.uniform(-b_2,b_2))
                    y_rand.append(random.uniform(-b_2,b_2))
                    angle_rand.append(random.uniform(float(0),float(360)))
                    xi_f.append(abs(random.normalvariate(xi_f_av,xi_f_stdev)))
                    delta_f.append(max(abs(random.normalvariate(delta_f_av,delta_f_stdev)),delta_f_min))
                    res_f.append(int(math.ceil(float(delta_f[fibrecount])/SplinePoints))) # fibre curvature resolution rounded to int

                    # debugging outputs
                    print(fibrecount) 
                    # print(x_rand[fibrecount])
                    # print(y_rand[fibrecount])
                    # print(angle_rand[fibrecount])

                    ## modelling
                    #part 1: sketch creation
                    m.ConstrainedSketch(name='__profile__', sheetSize=b_n)
                    a = m.rootAssembly

                    # create a single spline representing the sine curve
                    points = []
                    for x in range(int(-1.5*b_n), int(1.5*b_n), res_f[fibrecount]):
                        y = xi_f[fibrecount] * math.sin(2 * pi * x / delta_f[fibrecount])
                        points.append((x, y))

                    sketch = m.sketches['__profile__']
                    sketch.Spline(points=points)

                    # rotate & translate the spline
                    spline = sketch.geometry.values()[0]  # Assuming it's the only entity in the sketch
                    sketch.rotate(angle=angle_rand[fibrecount], centerPoint=(0.0, 0.0), objectList=[spline])
                    sketch.move(objectList=[spline], vector=(x_rand[fibrecount], y_rand[fibrecount]))
                    
                    # Define the coordinates of the rectangular area to keep (x_min, y_min) to (x_max, y_max)
                    x_min = -b_2
                    y_min = -b_2
                    x_max = b_2
                    y_max = b_2
                    
                    # find positions of the ends of the spline for removal purposes
                    pos1 = sketch.vertices[0].coords
                    pos2 = sketch.vertices[len(sketch.vertices)-1].coords
                    
                    # rectangle creation and setting to construction
                    sketch.rectangle(point1=(x_min,y_min),point2=(x_max,y_max))
                    sketch.setAsConstruction(objectList=(sketch.geometry.values()[1:5]))
                    
                    #trim curve to fit into RVE using the found endpoints
                    sketch.autoTrimCurve(curve1=sketch.geometry[2],point1=pos1)
                    sketch.autoTrimCurve(curve1=sketch.geometry[7],point1=pos2)

                    #finishing the sketch
                    m.Part(dimensionality=TWO_D_PLANAR, name=NoP, type=
                        DEFORMABLE_BODY)
                    m.parts[NoP].BaseWire(sketch=sketch)
                    del sketch
                    
                    # adding up the total length of fiber already "spun"
                    partlength = m.parts[NoP].edges[0].getSize()
                    AccLength += partlength*4
                    # creating sets 
                    pos1 = m.parts[NoP].vertices[0].pointOn
                    pos2 = m.parts[NoP].vertices[1].pointOn
                    pos3 = m.parts[NoP].edges[0].pointOn
                    
                    m.parts[NoP].Set(name='Start', vertices=
                        m.parts[NoP].vertices.findAt(pos1))    
                    m.parts[NoP].Set(name='End', vertices=
                        m.parts[NoP].vertices.findAt(pos2))    
                    m.parts[NoP].Set(edges=
                    m.parts[NoP].edges.findAt(pos3), name='fibre')
                    
                    # Assign section & beam orientation
                    m.parts[NoP].SectionAssignment(offset=0.0, 
                        offsetField='', offsetType=MIDDLE_SURFACE, region=
                        m.parts[NoP].sets['fibre'], sectionName=
                        'Circular_Beam', thicknessAssignment=FROM_SECTION)
                    m.parts[NoP].assignBeamSectionOrientation(method=
                        N1_COSINES, n1=(0.0, 0.0, -1.0), region=
                        m.parts[NoP].sets['fibre'])
                    
                    
                    #create mirroring planes
                    m.parts[NoP].DatumPlaneByPrincipalPlane(offset=0.0,
                        principalPlane=XZPLANE)
                    m.parts[NoP].DatumPlaneByPrincipalPlane(offset=0.0,
                        principalPlane=YZPLANE)
                    

                    # create mirrored copies
                    NoP_x = 'F'+str(fibrecount*4+1) # construct name of part -X
                    m.Part(name=NoP_x, objectToCopy=
                        m.parts[NoP])
                    m.parts[NoP_x].Mirror(keepOriginal=OFF, mirrorPlane=
                        m.parts[NoP_x].datums[5])
                    
                    NoP_y = 'F'+str(fibrecount*4+2) # construct name of part -Y
                    m.Part(name=NoP_y, objectToCopy=
                        m.parts[NoP])
                    m.parts[NoP_y].Mirror(keepOriginal=OFF, mirrorPlane=
                        m.parts[NoP_y].datums[6])
                    
                    NoP_xy = 'F'+str(fibrecount*4+3) # construct name of part -XY
                    m.Part(name=NoP_xy, objectToCopy=
                        m.parts[NoP_x])
                    m.parts[NoP_xy].Mirror(keepOriginal=OFF, mirrorPlane=
                        m.parts[NoP_xy].datums[6])

                    fibrecount += 1 # iterate partnumber
                    # sinusoidal part generation finished at this point --> from now on we work in assembly   

            # randomwalk fibre generation
            elif FibreType == 'RandomWalk':
                fibrecount = 0
                AccLength = 0
                ang_dev = []
                ang_dev.append(min(max(random.normalvariate(ang_av,ang_stdev),ang_min),ang_max))
                angle_rand = []
                # Define the coordinates of the rectangular area to keep (x_min, y_min) to (x_max, y_max)
                x_min = -float(b_n)/2
                y_min = -float(b_n)/2
                x_max = float(b_n)/2
                y_max = float(b_n)/2
                #defining the start point as the bottom left corner to hold still in PBCs
                points = []
                x = x_min
                y = y_min
                points.append((x,y))

                Puncture_LR = 0
                Puncture_UD = 0

                ang_n = 45 # start angle at 45°
                while AccLength<(TotalFibreLength): # create a bunch of fibers
                    
                    NoP = 'F'+str(fibrecount) # construct name of part
                    
                    
                    # sketch creation
                    m.ConstrainedSketch(name='__profile__', sheetSize=b_n)
                    if len(points)>1: # this makes sure the angle of a newly started fibre cannot shift outside the boundaries
                        ang_n = ang_n+random.uniform(-ang_dev[fibrecount],ang_dev[fibrecount]) #variate the angle between segments
                    x_old = x
                    y_old = y
                    dx = math.cos(ang_n)*l_step
                    dy = math.sin(ang_n)*l_step
                    x = x_old+math.cos(ang_n)*l_step # create next point x
                    y = y_old+math.sin(ang_n)*l_step*Alg # create next point y
                    new_ang_n = math.atan((y-y_old)/(x-x_old)) #account for the parameter Alg, calculate the made real angle
                    points.append((x,y))
                    angle_rand.append(math.degrees(new_ang_n)%360)    
                    
                    # check if new point penetrated any rve boundaries & if so, create a part out of existing splinepoints 
                    # & lay foundation for next part
                    if abs(x) > x_max or abs(y)> y_max:
                    #Figure out where it crossed exactly and trim the last point to be on the boundary
                        if x > x_max:
                            x_end = x_max
                            y_end = (x_max-x_old)/(x-x_old)*(y-y_old)+y_old
                            points[-1] = (x_end,y_end)
                            x_new = x_min # define startpoints for next fibre
                            y_new = y_end
                            Puncture_LR += 1
                        elif x < x_min:
                            x_end = x_min
                            y_end = (x_min-x_old)/(x-x_old)*(y-y_old)+y_old
                            points[-1] = (x_end,y_end)
                            x_new = x_max
                            y_new = y_end
                            Puncture_LR += 1
                        if y > y_max:
                            y_end = y_max
                            x_end = (y_max-y_old)/(y-y_old)*(x-x_old)+x_old
                            points[-1] = (x_end,y_end)
                            x_new = x_end
                            y_new = y_min
                            Puncture_UD += 1
                        elif y < y_min:
                            y_end = y_min
                            x_end = (y_min-y_old)/(y-y_old)*(x-x_old)+x_old
                            points[-1] = (x_end,y_end)
                            x_new = x_end
                            y_new = y_max
                            Puncture_UD += 1
                        x = x_new
                        y = y_new

                        print(AccLength)

                        sketch = m.sketches['__profile__']
                        sketch.Spline(points=points)         
                        points = [] # clear out old path 
                        points.append((x,y))           
                        m.Part(dimensionality=TWO_D_PLANAR, name=NoP, type=DEFORMABLE_BODY)
                        m.parts[NoP].BaseWire(sketch=sketch)
                        del sketch
                        # add up partlength
                        partlength = m.parts[NoP].edges[0].getSize()
                        AccLength += partlength
                        # creating sets 
                        pos1 = m.parts[NoP].vertices[0].pointOn
                        pos2 = m.parts[NoP].vertices[1].pointOn
                        pos3 = m.parts[NoP].edges[0].pointOn
                        
                        m.parts[NoP].Set(name='Start', vertices=
                            m.parts[NoP].vertices.findAt(pos1))    
                        m.parts[NoP].Set(name='End', vertices=
                            m.parts[NoP].vertices.findAt(pos2))    
                        m.parts[NoP].Set(edges=
                        m.parts[NoP].edges.findAt(pos3), name='fibre')

                        # Assign section & beam orientation
                        m.parts[NoP].SectionAssignment(offset=0.0,
                            offsetField='', offsetType=MIDDLE_SURFACE, region=
                            m.parts[NoP].sets['fibre'], sectionName=
                            'Circular_Beam', thicknessAssignment=FROM_SECTION)
                        m.parts[NoP].assignBeamSectionOrientation(method=
                            N1_COSINES, n1=(0.0, 0.0, -1.0), region=
                            m.parts[NoP].sets['fibre'])

                        ang_dev.append(min(max(random.normalvariate(ang_av,ang_stdev),ang_min),ang_max))
                        fibrecount +=1
            FLS = []
            for p in range(len(m.parts)): 
                FLS.append('F'+str(p))
            random.shuffle(FLS)

            mdb.customData.SaveFLS = FLS #save randomized fibre list if crosslinks want to be regenerated in the same way
            mdb.customData.savefloats.append(AccLength) #save actual generated fibrelength

            # create angle histogram
            ang_hist = [0,0,0,0,0,0,0,0]
            for a in angle_rand:
                b = a%180
                if (0<b<=22.5):
                    ang_hist[4]+=1
                elif (22.5<b<=45):
                    ang_hist[5]+=1
                elif (45<b<=67.5):
                    ang_hist[6]+=1
                elif (67.5<b<=90):
                    ang_hist[7]+=1
                elif (90<b<=112.5):
                    ang_hist[0]+=1
                elif (112.5<b<=135):
                    ang_hist[1]+=1
                elif (135<b<=157.5):
                    ang_hist[2]+=1
                elif (157.5<b<=180):
                    ang_hist[3]+=1

            # account for mirrored copies in sinuosidal fibres
            if FibreType == 'Sinusoidal':
                ang_hist_new = [0,0,0,0,0,0,0,0]
                for i in range(len(ang_hist)):
                    ang_hist_new[i] = (ang_hist[i]+ang_hist[len(ang_hist)-1-i])*2
                ang_hist = ang_hist_new

            mdb.customData.saveLists = [ang_hist]
            mean_ang_hist = sum(ang_hist) / len(ang_hist) # Calculates the mean of the values
            print(mean_ang_hist) 
            squared_deviations = [(x - mean_ang_hist) ** 2 for x in ang_hist] 
            stdev_ang_hist = math.sqrt(sum(squared_deviations) / (len(ang_hist)-1))
            print(stdev_ang_hist) # Calculates the standard deviation of the values 
            v_iso = 1 - 1/mean_ang_hist * stdev_ang_hist  # Calculate the isotropy coefficient
            print('Angle Histogram:'+str(ang_hist))

            
            
            
            
            
            
            ## start of crosslink & PBC generation

            a = m.rootAssembly
            a.DatumCsysByDefault(CARTESIAN)

            # creating reference points for deformation loading
            a.ReferencePoint(point=(-b_2, -b_2, 0.0))
            a.ReferencePoint(point=(b_2, -b_2, 0.0))
            a.ReferencePoint(point=(-b_2, b_2, 0.0))
            a.features.changeKey(fromName='RP-1', toName= 'Origin')
            a.features.changeKey(fromName='RP-2', toName= 'X-Dir')
            a.features.changeKey(fromName='RP-3', toName= 'Y-Dir')
            a.Set(name='Origin', referencePoints=(a.referencePoints[2], ))
            a.Set(name='X-Dir', referencePoints=(a.referencePoints[3], ))
            a.Set(name='Y-Dir', referencePoints=(a.referencePoints[4], ))

            FL = []
            FLS = mdb.customData.SaveFLS
            for p in range(len(m.parts)):
                FL.append('F'+str(p))

                    # Partition fibers in half until edglength is lower than res_fp, always partition them at least once
                edgeLength = m.parts[FL[p]].edges[0].getSize(printResults=FALSE)
                m.parts[FL[p]].PartitionEdgeByParam(edges=
                    m.parts[FL[p]].edges, parameter=0.5)
                NodeCounter = 2
                while edgeLength > res_fp:
                    m.parts[FL[p]].PartitionEdgeByParam(edges=
                        m.parts[FL[p]].edges, parameter=0.5)
                    edgeLength = m.parts[FL[p]].edges[0].getSize(printResults=FALSE)
                    NodeCounter *=2

                    # Set element type and mesh the part
                m.parts[FL[p]].setElementType(elemTypes=(ElemType(
                            elemCode=B22, elemLibrary=STANDARD), ), regions=(
                                m.parts[FL[p]].sets['fibre']))
                m.parts[FL[p]].seedPart(deviationFactor=0.1,
                            minSizeFactor=0.1, size=res_fp)
                m.parts[FL[p]].generateMesh()

                    # Create nodal set without start & end node inside (for crosslinks)
                if NodeCounter == 2:
                    m.parts[FL[p]].SetFromNodeLabels(name = 'InternalNodes',nodeLabels=[2])
                else:
                    m.parts[FL[p]].SetFromNodeLabels(name = 'InternalNodes',nodeLabels=(
                        tuple(range(2,NodeCounter))))
                
                    # Load part into Assembly Instance
                a.Instance(dependent=ON, name=FL[p], part= m.parts[FL[p]])

            
                
            # calculate new interaction range Int_t based on a consistent mean segment length l_s
            #interaction range --> how many neighboring fibers get taken into account
            if l_s_mean: # if l_s_mean is not equal to False, it is executed
                print('Consistent mean segment length activated l_s_mean:'+str(l_s_mean))
                t_I = d_f**2/(8*(1-por))*math.pi**2/l_s_mean # interaction thickness of one fiber. Eqn. 10 from DOI: 10.1039/c7sm01241a
                t_F = float(t_n)/float((fibrecount+1)) #smeared thickness of one fiber
                Int_t = int(round(t_I/t_F))
                l_s_mean_actual = d_f**2/(8*(1-por))*math.pi**2/(Int_t*t_F) #calculating the actual l_s_mean
                print('Calculated interaction range for crosslink generation:'+ str(Int_t))
                print('Calculated actual mean segment length:'+ str(l_s_mean_actual))
                print('Interaction thickness: '+str(t_I))



            # creating the crosslinks in shuffled stack
            update_threshold = 0.0
            Comptot = len(a.instances)-Int_t
            for f in range(Int_t,len(a.instances)):
                F_slave = FLS[f]
                perc_done = (f-Int_t)/float(Comptot)
                if perc_done >= update_threshold:
                    print(str(int(perc_done*100))+'%')
                    update_threshold = update_threshold+0.1
                F_master = FLS[f-Int_t]
                a.SetByBoolean(name='Temp', sets=(a.sets[F_master+'.InternalNodes'],))
                for g in range(f-Int_t+1,f):
                    F_master = FLS[g]
                    a.SetByBoolean(name='Temp', sets=(
                        a.sets['Temp'],a.sets[F_master+'.InternalNodes'],))
                
                NoT0 = FLS[f-Int_t]+'---'+FLS[f-1]
                a.SetByBoolean(name=NoT0, sets=(a.sets['Temp'],))
                del a.sets['Temp']
                NoT1 = F_slave+'-'+NoT0
                m.Tie(adjust=OFF, main=
                    a.sets[NoT0], name=NoT1, positionTolerance=float(res_fp*0.71), 
                    positionToleranceMethod=SPECIFIED, secondary=
                    a.sets[F_slave+'.InternalNodes'], thickness=OFF, tieRotations=ON)

            # creating the Periodic Boundary Conditions
            S_E = ['.Start','.End']
            S_E_s = ['S','E']
            if FibreType == 'Sinusoidal':
                C_O = 0
                for k in range(0,len(FL),4):
                    
                    for l in range(len(S_E)):
                        Pos = a.sets[FL[k]+S_E[l]].vertices[0].pointOn[0]
                        CoeffLR = Pos[0]/(b_2)
                        CoeffUD = Pos[1]/(b_2)
                        # Left-Right linking
                        if (abs(CoeffLR)==1):
                            m.Equation(name='PBC'+str(k)+S_E_s[l]+'-LR-X', terms=(
                                (CoeffLR, FL[k]+S_E[l], 1), (-CoeffLR, FL[k+2]+S_E[l], 1),
                                (-1.0, 'X-Dir', 1), (1.0, 'Origin', 1)))
                            m.Equation(name='PBC'+str(k)+S_E_s[l]+'-LR-Y', terms=(
                                (CoeffLR, FL[k]+S_E[l], 2), (-CoeffLR, FL[k+2]+S_E[l], 2),
                                (-1.0, 'X-Dir', 2), (1.0, 'Origin', 2)))
                            m.Equation(name='PBC'+str(k)+S_E_s[l]+'-Rot', terms=(
                                (CoeffLR, FL[k]+S_E[l], 6), (-CoeffLR, FL[k+2]+S_E[l], 6),))
                            
                            m.Equation(name='PBC'+str(k+1)+S_E_s[l]+'-LR-X', terms=(
                                (CoeffLR, FL[k+1]+S_E[l], 1), (-CoeffLR, FL[k+3]+S_E[l], 1),
                                (-1.0, 'X-Dir', 1), (1.0, 'Origin', 1)))
                            m.Equation(name='PBC'+str(k+1)+S_E_s[l]+'-LR-Y', terms=(
                                (CoeffLR, FL[k+1]+S_E[l], 2), (-CoeffLR, FL[k+3]+S_E[l], 2),
                                (-1.0, 'X-Dir', 2), (1.0, 'Origin', 2)))
                            m.Equation(name='PBC'+str(k+1)+S_E_s[l]+'-Rot', terms=(
                                (CoeffLR, FL[k+1]+S_E[l], 6), (-CoeffLR, FL[k+3]+S_E[l], 6),))
                            
                            # find the LR-set closest to origin
                            if abs(CoeffUD) > C_O:
                                C_O = abs(CoeffUD)
                                if (CoeffLR == 1 and CoeffUD >0):
                                    C_set = FL[k+3]+S_E[l]
                                elif (CoeffLR == 1 and CoeffUD <0):
                                    C_set = FL[k+2]+S_E[l]
                                elif (CoeffLR == -1 and CoeffUD >0):
                                    C_set = FL[k+1]+S_E[l]
                                elif (CoeffLR == -1 and CoeffUD <0):
                                    C_set = FL[k]+S_E[l]

                        # Up-Down linking
                        elif (abs(CoeffUD)==1):
                            m.Equation(name='PBC'+str(k)+S_E_s[l]+'-UD-X', terms=(
                                (CoeffUD, FL[k]+S_E[l], 1), (-CoeffUD, FL[k+1]+S_E[l], 1),
                                (-1.0, 'Y-Dir', 1), (1.0, 'Origin', 1)))
                            m.Equation(name='PBC'+str(k)+S_E_s[l]+'-UD-Y', terms=(
                                (CoeffUD, FL[k]+S_E[l], 2), (-CoeffUD, FL[k+1]+S_E[l], 2),
                                (-1.0, 'Y-Dir', 2), (1.0, 'Origin', 2)))
                            m.Equation(name='PBC'+str(k)+S_E_s[l]+'-Rot', terms=(
                                (CoeffUD, FL[k]+S_E[l], 6), (-CoeffUD, FL[k+1]+S_E[l], 6),))
                            
                            m.Equation(name='PBC'+str(k+1)+S_E_s[l]+'-UD-X', terms=(
                                (CoeffUD, FL[k+2]+S_E[l], 1), (-CoeffUD, FL[k+3]+S_E[l], 1),
                                (-1.0, 'Y-Dir', 1), (1.0, 'Origin', 1)))
                            m.Equation(name='PBC'+str(k+1)+S_E_s[l]+'-UD-Y', terms=(
                                (CoeffUD, FL[k+2]+S_E[l], 2), (-CoeffUD, FL[k+3]+S_E[l], 2),
                                (-1.0, 'Y-Dir', 2), (1.0, 'Origin', 2)))
                            m.Equation(name='PBC'+str(k+1)+S_E_s[l]+'-Rot', terms=(
                                (CoeffUD, FL[k+2]+S_E[l], 6), (-CoeffUD, FL[k+3]+S_E[l], 6),))
                            

            elif FibreType == 'RandomWalk':
                C_set = 'F0.Start'
                for k in range(0,len(FL)-1):
                    Pos = a.sets[FL[k]+'.End'].vertices[0].pointOn[0]
                    CoeffLR = Pos[0]/(b_2)
                    CoeffUD = Pos[1]/(b_2)
                    # Left-Right linking
                    if (abs(CoeffLR)==1):
                        m.Equation(name='PBC'+str(k)+str(k+1)+'-LR-X', terms=(
                            (CoeffLR, FL[k]+'.End', 1), (-CoeffLR, FL[k+1]+'.Start', 1),(
                                -1.0, 'X-Dir', 1), (1.0, 'Origin', 1)))
                        m.Equation(name='PBC'+str(k)+str(k+1)+'-LR-Y', terms=(
                            (CoeffLR, FL[k]+'.End', 2), (-CoeffLR, FL[k+1]+'.Start', 2),
                            (-1.0, 'X-Dir', 2), (1.0, 'Origin', 2)))
                        
                    elif (abs(CoeffUD)==1):
                        m.Equation(name='PBC'+str(k)+str(k+1)+'-UD-X', terms=(
                            (CoeffUD, FL[k]+'.End', 1), (-CoeffUD, FL[k+1]+'.Start', 1),
                            (-1.0, 'Y-Dir', 1), (1.0, 'Origin', 1)))
                        m.Equation(name='PBC'+str(k)+str(k+1)+'-UD-Y', terms=(
                            (CoeffUD, FL[k]+'.End', 2), (-CoeffUD, FL[k+1]+'.Start', 2),
                            (-1.0, 'Y-Dir', 2), (1.0, 'Origin', 2)))
                        
                    m.Equation(name='PBC'+str(k)+str(k+1)+'-Rot', terms=(
                        (CoeffLR, FL[k]+'.End', 6), (-CoeffLR, FL[k+1]+'.Start', 6),))



                        
            # creating the loadcases
            # initial holding
            m.DisplacementBC(amplitude=UNSET, createStepName='Initial', 
                distributionType=UNIFORM, fieldName='', localCsys=None, name='Origin', 
                region=a.sets['Origin'], u1=SET, u2=SET, u3=UNSET, 
                ur1=UNSET, ur2=UNSET, ur3=UNSET)
            m.DisplacementBC(amplitude=UNSET, createStepName='Initial', 
                distributionType=UNIFORM, fieldName='', localCsys=None, name='Y-Dir', 
                region=a.sets['Y-Dir'], u1=SET, u2=UNSET, u3=UNSET, 
                ur1=UNSET, ur2=UNSET, ur3=UNSET)
            m.DisplacementBC(amplitude=UNSET, createStepName='Initial', 
                distributionType=UNIFORM, fieldName='', localCsys=None, name='X-Dir', 
                region=a.sets['X-Dir'], u1=UNSET, u2=SET, 
                u3=UNSET, ur1=UNSET, ur2=UNSET, ur3=UNSET)

            # fix RBM by fixing LR set closest to origin (in RW, it is the start of the first fibre)
            m.DisplacementBC(amplitude=UNSET, createStepName='Initial', 
                distributionType=UNIFORM, fieldName='', localCsys=None, name='FixRBM', 
                region=a.sets[C_set], u1=SET, u2=SET, ur3=UNSET)

            # create Loadcases
            if FibreType == 'Sinusoidal':
                JobName_raw = date+'_S_'+str(b_n)+'_'+str(t_n)
            elif FibreType == 'RandomWalk':
                JobName_raw = date+'_R_'+str(b_n)+'x'+str(t_n)+'xlsm'+str(l_s_mean)+'_plast'+'_viso'+str(round(v_iso, 2))
            JobName_global = JobName_raw.replace('.','_')

            # Uniaxial
            if (LoadCase=='Uniaxial'):
                m.boundaryConditions['X-Dir'].setValuesInStep(stepName=
                    'LoadingStep-1', u1=eps11*float(b_n))
                # create Job
                # mdb.saveAs(JobName_global+'.cae')
                JobName = JobName_global+'_UniAx'
                with open(JobName+'_'+str(i)+'_fiberstackorder.txt', 'w') as file: # save list of fiber stack order
                    for item in FLS:
                        file.write("'"+item+"', ")
                    file.write('\ninteraction range: '+ str(Int_t))
                Desc = ('width = '+str(b_n)+', depth = '+str(t_n)+', mean segment length = '+str(l_s_mean)+
                        ', Resolution = '+str(res_fp)+', Fibrelength = '+str(AccLength)+', Number fibres= '+str(fibrecount)+', Plasticity = YES'+
                        ', Interaction thickness= '+str(t_I)+' ,isotropy coefficient= '+str(round(v_iso, 2))+', Alignment value= '+str(Alg))
                mdb.Job(atTime=None, contactPrint=OFF, description=Desc, echoPrint=OFF, 
                    explicitPrecision=SINGLE, getMemoryFromAnalysis=True, historyPrint=OFF, 
                    memory=90, memoryUnits=PERCENTAGE, model='Model-1', modelPrint=OFF, 
                    multiprocessingMode=DEFAULT, name=JobName+'_'+str(i), nodalOutputPrecision=SINGLE, 
                    numCpus=4, numDomains=4, numGPUs=0, queue=None, resultsFormat=ODB, scratch=
                    '', type=ANALYSIS, userSubroutine='', waitHours=0, waitMinutes=0)
                mdb.saveAs(JobName+'_'+str(i))
                mdb.jobs[JobName+'_'+str(i)].writeInput(consistencyChecking=OFF)
            # Planar Tension
            if (LoadCase=='PlanarTension'):
                m.boundaryConditions['Y-Dir'].setValuesInStep(stepName=
                    'LoadingStep-1', u2=SET)
                m.boundaryConditions['X-Dir'].setValuesInStep(stepName=
                    'LoadingStep-1', u1=eps11*float(b_n))
                # create Job
                # mdb.saveAs(JobName_global+'.cae')
                JobName = JobName_global+'_PlanTen'
                with open(JobName+'_'+str(i)+'_fiberstackorder.txt', 'w') as file: # save list of fiber stack order
                    for item in FLS:
                        file.write("'"+item+"', ")
                    file.write('\ninteraction range: '+ str(Int_t))
                Desc = ('width = '+str(b_n)+', depth = '+str(t_n)+', mean segment length = '+str(l_s_mean)+
                        ', Resolution = '+str(res_fp)+', Fibrelength = '+str(AccLength)+', Number fibres= '+str(fibrecount)+', Plasticity = YES'+
                        ', Interaction thickness= '+str(t_I)+' ,isotropy coefficient= '+str(round(v_iso, 2))+', Alignment value= '+str(Alg))
                mdb.Job(atTime=None, contactPrint=OFF, description=Desc, echoPrint=OFF, 
                    explicitPrecision=SINGLE, getMemoryFromAnalysis=True, historyPrint=OFF, 
                    memory=90, memoryUnits=PERCENTAGE, model='Model-1', modelPrint=OFF, 
                    multiprocessingMode=DEFAULT, name=JobName+'_'+str(i), nodalOutputPrecision=SINGLE, 
                    numCpus=4, numDomains=4, numGPUs=0, queue=None, resultsFormat=ODB, scratch=
                    '', type=ANALYSIS, userSubroutine='', waitHours=0, waitMinutes=0)
                mdb.saveAs(JobName+'_'+str(i))
                mdb.jobs[JobName+'_'+str(i)].writeInput(consistencyChecking=OFF)

            # Simple Shear
            if (LoadCase=='SimpleShear'):
                m.boundaryConditions['Y-Dir'].setValuesInStep(stepName=
                    'LoadingStep-1', u1=eps12*2*float(b_n))
                m.boundaryConditions['X-Dir'].setValuesInStep(stepName=
                    'LoadingStep-1', u1=SET)
                # create Job
                # mdb.saveAs(JobName_global+'.cae')
                JobName = JobName_global+'_SimpleShear'
                with open(JobName+'_'+str(i)+'_fiberstackorder.txt', 'w') as file: # save list of fiber stack order
                    for item in FLS:
                        file.write("'"+item+"', ")
                    file.write('\ninteraction range: '+ str(Int_t))
                Desc = ('width = '+str(b_n)+', depth = '+str(t_n)+', mean segment length = '+str(l_s_mean)+
                        ', Resolution = '+str(res_fp)+', Fibrelength = '+str(AccLength)+', Number fibres= '+str(fibrecount)+', Plasticity = YES'+
                        ', Interaction thickness= '+str(t_I)+' ,isotropy coefficient= '+str(round(v_iso, 2))+', Alignment value= '+str(Alg))
                mdb.Job(atTime=None, contactPrint=OFF, description=Desc, echoPrint=OFF, 
                    explicitPrecision=SINGLE, getMemoryFromAnalysis=True, historyPrint=OFF, 
                    memory=90, memoryUnits=PERCENTAGE, model='Model-1', modelPrint=OFF, 
                    multiprocessingMode=DEFAULT, name=JobName+'_'+str(i), nodalOutputPrecision=SINGLE, 
                    numCpus=4, numDomains=4, numGPUs=0, queue=None, resultsFormat=ODB, scratch=
                    '', type=ANALYSIS, userSubroutine='', waitHours=0, waitMinutes=0)
                mdb.saveAs(JobName+'_'+str(i))
                mdb.jobs[JobName+'_'+str(i)].writeInput(consistencyChecking=OFF)
    print('Angle Histogram:'+str(ang_hist))
    print('isotropy coefficient'+str(v_iso))