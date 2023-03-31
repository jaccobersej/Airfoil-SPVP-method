import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import math
from matplotlib import path


#Kendte Variable:
angles = [0, 2.5, 5, 7.5, 10]

#Vindhastighed
Vinf = 10.8                    
for i in angles:
    #Angle of attack (Vinkel på vingen i forhold til vinden)
    AoA = i                     


    #Omregn vinkel til radianer:
    alpha = AoA*(np.pi/180)






    #LAV AIRFOIL


    #Load airfoil koordinater lavet i loggerPro ved billed analyse
    airfoilCoords = np.array(pd.read_excel("Airfoil140.xlsx")).T
    xBoundary = (airfoilCoords[0]-min(airfoilCoords[0]))*4.69071171
    yBoundary = (airfoilCoords[1]+(min(airfoilCoords[1])*-1))*4.69071171
    boundaryNumber = len(xBoundary)
    panelNumber = boundaryNumber - 1


    #Check om clockwise:
    edge = np.zeros(panelNumber)


    for i in range(panelNumber):
        edge[i] = (xBoundary[i+1]-xBoundary[i])*(yBoundary[i+1]+yBoundary[i])

    sumEdge = np.sum(edge)  


    #Hvis sumEdge < 0, må den være rækkefølgen være counterclockwise og skal derfor flippes
    if (sumEdge < 0):
        xBoundary = np.flipud(xBoundary)
        yBoundary = np.flipud(yBoundary)
       
    #initialize tomme arrays (Geometric quants)
    xControlPoints = np.zeros(panelNumber)
    yControlPoints = np.zeros(panelNumber)
    panelLengths = np.zeros(panelNumber)
    angleBetweenXandPanel = np.zeros(panelNumber)


    #Find geometric quants
    for i in range(panelNumber):
        xControlPoints[i] = 0.5*(xBoundary[i]+xBoundary[i+1])
        yControlPoints[i] = 0.5*(yBoundary[i]+yBoundary[i+1])
        dx = xBoundary[i+1] - xBoundary[i]
        dy = yBoundary[i+1] - yBoundary[i]
        panelLengths[i] = (dx**2+dy**2)**0.5
        angleBetweenXandPanel[i] = math.atan2(dy, dx)
        if (angleBetweenXandPanel[i] < 0):
            angleBetweenXandPanel[i] = angleBetweenXandPanel[i] + 2*np.pi
       
    #relevante vinkler
    delta = angleBetweenXandPanel + (np.pi/2)
    beta = delta - alpha
    beta[beta > 2*np.pi] = beta[beta > 2*np.pi] - 2*np.pi


    #FIND PANEL STRENGTHS


    def SourceGeometricIntegrals(xControl, yControl, xBounds, yBounds, phi, pLengths):
       
        #Init arrays
        I_values = np.zeros([panelNumber, panelNumber])
        J_values = np.zeros([panelNumber, panelNumber])
       
        #Beregn integraler
        for i in range(panelNumber):
            for j in range(panelNumber):
                if (j != i):
                    #Delværdier, som bruges til den samlede udregning. Skrevet op for sig selv for at gøre det mere overskueligt
                    a = -(xControl[i]-xBounds[j])*np.cos(phi[j])-(yControl[i]-yBounds[j])*np.sin(phi[j])
                    b = (xControl[i]-xBounds[j])**2 + (yControl[i]-yBounds[j])**2
                    cNormal = np.sin(phi[i]-phi[j])
                    dNormal = -(xControl[i]-xBounds[j])*np.sin(phi[i])+(yControl[i]-yBounds[j])*np.cos(phi[i])
                    cTangent = -np.cos(phi[i]-phi[j])
                    dTangent = (xControl[i]-xBounds[j])*np.cos(phi[i])+(yControl[i]-yBounds[j])*np.sin(phi[i])
                    e_value = np.sqrt(b-a**2)
                    if (e_value == 0 or np.iscomplex(e_value) or np.isnan(e_value) or np.isinf(e_value)):
                        I_values[i,j] = 0
                        J_values[i,j] = 0
                    else:
                        #Beregn I
                        term1 = 0.5 * cNormal * np.log((pLengths[j]**2 + 2*a*pLengths[j] + b)/b)
                        term2 = ((dNormal - a*cNormal)/e_value) * (math.atan2((pLengths[j]+ a), e_value)-math.atan2(a, e_value))
                        I_values[i,j] = term1 + term2
                       
                        #Beregn J
                        term1 = 0.5 * cTangent * np.log((pLengths[j]**2 + 2*a*pLengths[j] + b)/b)
                        term2 = ((dTangent - a*cTangent)/e_value) * (math.atan2((pLengths[j]+ a), e_value)-math.atan2(a, e_value))
                        J_values[i,j] = term1 + term2
       
        return I_values, J_values


    #geometric integrals K og L
    def VortexGeometricIntegrals(xControl, yControl, xBounds, yBounds, phi, pLengths):
       
       
       
        #Initialize array:
        k_values = np.zeros([panelNumber, panelNumber])
        L_values = np.zeros([panelNumber, panelNumber])
       
        #Beregn integrale
        for i in range(panelNumber):
            for j in range(panelNumber):
                if (j != i):
                    #Delværdier, som bruges til den samlede udregning. Skrevet op for sig selv for at gøre det mere overskueligt
                    a = -(xControl[i]-xBounds[j])*np.cos(phi[j])-(yControl[i]-yBounds[j])*np.sin(phi[j])
                    b = (xControl[i]-xBounds[j])**2 + (yControl[i]-yBounds[j])**2
                    cNormal = -np.cos(phi[i]-phi[j])
                    dNormal = (xControl[i]-xBounds[j])*np.cos(phi[i])+(yControl[i]-yBounds[j])*np.sin(phi[i])
                    cTangent = np.sin(phi[j]-phi[i])
                    dTangent = (xControl[i]-xBounds[j])*np.sin(phi[i])-(yControl[i]-yBounds[j])*np.cos(phi[i])
                    e_value = np.sqrt(b-(a**2))
                    if (e_value == 0 or np.iscomplex(e_value) or np.isnan(e_value) or np.isinf(e_value)):
                        k_values[i,j] = 0
                        L_values[i,j] = 0
                    else:
                        #Beregn k
                        term1 = 0.5 * cNormal * np.log((pLengths[j]**2 + 2*a*pLengths[j] + b)/b)
                        term2 = ((dNormal - a*cNormal)/e_value) * (math.atan2((pLengths[j]+ a), e_value)-math.atan2(a, e_value))
                        k_values[i,j] = term1 + term2
                       
                        #Beregn L
                        term1 = 0.5 * cTangent * np.log((pLengths[j]**2 + 2*a*pLengths[j] + b)/b)
                        term2 = ((dTangent - a*cTangent)/e_value) * (math.atan2((pLengths[j]+ a), e_value)-math.atan2(a, e_value))
                        L_values[i,j] = term1 + term2
                       
                    if (np.iscomplex(k_values[i,j]) or np.isnan(k_values[i,j]) or np.isinf(k_values[i,j])):      
                        k_values[i,j] = 0    
                           
        return k_values, L_values


    I, J = SourceGeometricIntegrals(xControlPoints, yControlPoints, xBoundary, yBoundary, angleBetweenXandPanel, panelLengths)
    K, L = VortexGeometricIntegrals(xControlPoints, yControlPoints, xBoundary, yBoundary, angleBetweenXandPanel, panelLengths)


    #Lav matrix A
    A = np.zeros([panelNumber, panelNumber])
    for i in range(panelNumber):
        for j in range(panelNumber):
            if i == j:
                A[i,j] = np.pi
            else:
                A[i,j] = I[i,j]    
               
    #Ekstra kolonne til A matrix:
    newACVector = np.zeros((panelNumber, 1))
    A = np.hstack((A, newACVector))
    for i in range(panelNumber):
        A[i, panelNumber] = -sum(K[i,:])


    #Forstør matrix til kutta condition:
    newARVector = np.zeros((1, panelNumber+1))
    A = np.vstack((A, newARVector))
    for j in range(panelNumber):
        A[panelNumber, j] = J[0, j] + J[panelNumber-1, j]
    A[panelNumber, panelNumber] = -(sum(L[0,:] + L[panelNumber-1,:])) + 2*np.pi


    #Lav matrix B
    b = np.zeros(panelNumber)
    for i in range(panelNumber):
        b[i] = -Vinf*2*np.pi*np.cos(beta[i])


    b = np.append(b,-Vinf*2*np.pi*(np.sin(beta[0]) + np.sin(beta[panelNumber-1])))


    #Compute vortex strengths via linalg equation system solve
    result = np.linalg.solve(A,b)

    #Unpack source strengths og vortex strength for resultatet
    lambdaCollection = result[0:len(result)-1]
    gamma = result[len(result)-1]
       
    #Beregn den tangentiale hastighed og trykkoefficienter 
    Vt = np.zeros(panelNumber)                                                          
    Cp = np.zeros(panelNumber)                                                          
    for i in range(panelNumber):                                                        
        term1 = Vinf*np.sin(beta[i])
        term2 = (1/(2*np.pi))*sum(lambdaCollection*J[i,:])                        
        term3 = gamma/2
        term4 = -(gamma/(2*np.pi))*sum(L[i,:])
       
       
        Vt[i] = term1 + term2 + term3+ term4                        
        Cp[i] = 1 - (Vt[i]/Vinf)**2  
       
       
    #Bryd trykkoefficienterne op i deres normal og axiale kraft koefficienter
    CN = -Cp*panelLengths*np.sin(beta)                                                        
    CA = -Cp*panelLengths*np.cos(beta)        
    print(str(AoA) + "°")
    #Beregne lift koefficient via ovenstående
    CL = sum(CN*np.cos(alpha)) - sum(CA*np.sin(alpha))      
    print("Lift koefficient fundet via tryk koefficienter:", CL)


    #COMPUTE LIFT COEFFICIENT


        #Compute via Kutta-Joukowski:
    liftCoeff = (2 * sum(gamma*panelLengths))/(Vinf)


    print("Lift koefficient fundet via Kutta-Joukowski udledning:", liftCoeff)
        #Print
       
       


    #Beregn værdier til streamlines
    """
    def streamlineSPM(xP, yP, xB, yB, phi, pLengths):
       
        mx = np.zeros(panelNumber)
        my = np.zeros(panelNumber)
       
        #beregn integrale:
        for j in range(panelNumber):
            a = -(xP-xB[j])*np.cos(phi[j]) - (yP-yB[j])*np.sin(phi[j])
            b = (xP-xB[j])**2 + (yP-yB[j])**2
            cX = np.sin(phi[j])
            dX = -(yP - yB[j])
            cY = -np.cos(phi[j])
            dY = xP - xB[j]
            e_value = np.sqrt(b-a**2)
            if (e_value == 0 or np.iscomplex(e_value) or np.isnan(e_value) or np.isinf(e_value)):
                mx[j] = 0
                my[j] = 0
            else:
                #Beregn k
                term1 = 0.5 * cX * np.log((pLengths[j]**2 + 2*a*pLengths[j] + b)/b)
                term2 = ((dX - a*cX)/e_value) * (math.atan2((pLengths[j]+ a), e_value)-math.atan2(a, e_value))
                mx[j] = term1 + term2
               
                #Beregn L
                term1 = 0.5 * cY * np.log((pLengths[j]**2 + 2*a*pLengths[j] + b)/b)
                term2 = ((dY - a*cY)/e_value) * (math.atan2((pLengths[j]+ a), e_value)-math.atan2(a, e_value))
                my[j] = term1 + term2
       
        return mx, my


    def streamlineVPM(xP, yP, xB, yB, phi, pLengths):
       
        nx = np.zeros(panelNumber)
        ny = np.zeros(panelNumber)
       
        #beregn integrale:
        for j in range(panelNumber):
            a = -(xP-xB[j])*np.cos(phi[j]) - (yP-yB[j])*np.sin(phi[j])
            b = (xP-xB[j])**2 + (yP-yB[j])**2
            cX = -np.cos(phi[j])
            dX = xP - xB[j]
            cY = -np.sin(phi[j])
            dY = yP - yB[j]
            e_value = np.sqrt(b-a**2)
            if (e_value == 0 or np.iscomplex(e_value) or np.isnan(e_value) or np.isinf(e_value)):
                nx[j] = 0
                ny[j] = 0
            else:
                #Beregn k
                term1 = 0.5 * cX * np.log((pLengths[j]**2 + 2*a*pLengths[j] + b)/b)
                term2 = ((dX - a*cX)/e_value) * (math.atan2((pLengths[j]+ a), e_value)-math.atan2(a, e_value))
                nx[j] = term1 + term2
               
                #Beregn L
                term1 = 0.5 * cY * np.log((pLengths[j]**2 + 2*a*pLengths[j] + b)/b)
                term2 = ((dY - a*cY)/e_value) * (math.atan2((pLengths[j]+ a), e_value)-math.atan2(a, e_value))
                ny[j] = term1 + term2
       
        return nx, ny


    #Grid parametre
    nGridX = 500                                                              
    nGridY = 500                                                      
    xVals  = [min(xBoundary)-0.5, max(xBoundary)+0.5]                                    
    yVals  = [min(yBoundary)-0.6, max(yBoundary)+0.6]                                


    #Streamline parametre
    slPct  = 20                                                                
    Ysl    = np.linspace(yVals[0]+0.01,yVals[1]-0.01,int((slPct/100)*nGridY))            
    Xsl    = xVals[0]*np.ones(len(Ysl))                                        
    XYsl   = np.vstack((Xsl.T,Ysl.T)).T                                        


    print(XYsl.shape)
    #Lav grid points
    Xgrid  = np.linspace(xVals[0],xVals[1],nGridX)                            
    Ygrid  = np.linspace(yVals[0],yVals[1],nGridY)
    #xValuesForStack =
    #stack = np.hstack((Xgrid[1:len(Xgrid)-1].reshape(98,1), Ygrid[1:len(Ygrid)-1].reshape(98,1)))    
    #print(stack.shape)            
    XX, YY = np.meshgrid(Xgrid, Ygrid)    
       
    #Init arrays
    Vx = np.zeros([nGridX, nGridY])
    Vy = np.zeros([nGridX, nGridY])


    #Path for at checke om punkter er inden i vingen:
    AF = np.vstack((xBoundary.T, yBoundary.T)).T
    afPath = path.Path(AF)


    for m in range(nGridX):
        print("m: %i" % m)
        for n in range(nGridY):
            xPoint = XX[m, n]
            yPoint = YY[m, n]
            mx, my = streamlineSPM(xPoint, yPoint, xBoundary, yBoundary, angleBetweenXandPanel, panelLengths)
            nx, ny = streamlineVPM(xPoint, yPoint, xBoundary, yBoundary, angleBetweenXandPanel, panelLengths)
           
           
            if afPath.contains_points([(xPoint, yPoint)]):
                Vx[m,n] = 0
                Vy[m,n] = 0
            else:
                Vx[m,n] = (Vinf*np.cos(alpha)+ sum(lambdaCollection*mx/(2*np.pi))+ sum(-gamma*nx/(2*np.pi)))
                Vy[m,n] = (Vinf*np.sin(alpha)+ sum(lambdaCollection*my/(2*np.pi))+ sum(-gamma*ny/(2*np.pi)))
               
               
    Vxy = np.sqrt(Vx**2 + Vy**2)
    pressureCoeffXY = 1 - (Vxy/Vinf)**2
       
        """
    #plot af tryk koefficient normaler
    fig = plt.figure(1)                                                
    plt.cla()                                          
    Cps = np.absolute(Cp*0.15)                                                
    X = np.zeros(2)                                                        
    Y = np.zeros(2)                                                    
    for i in range(len(Cps)):                                        
        #print(Cp[i])
        X[0] = xControlPoints[i]                                                          
        X[1] = xControlPoints[i] + Cps[i]*np.cos(delta[i])                              
        Y[0] = yControlPoints[i]                                                            
        Y[1] = yControlPoints[i] + Cps[i]*np.sin(delta[i])                                  
       
        if (Cp[i] < 0):                                                      
            plt.plot(X,Y,'r-')                                                  
        elif (Cp[i] >= 0):                                                    
            plt.plot(X,Y,'b-')                                                  
    plt.fill(xBoundary,yBoundary,'k')                                                    
    plt.xlabel('X Units')                                                    
    plt.ylabel('Y Units')                                                      
    plt.gca().set_aspect('equal')                                            
    plt.show()                  
    """

    #Streamlines plot
    fig = plt.figure(2)                                                         # Create figure
    plt.cla()                                                                   # Get ready for plotting
    np.seterr(under="ignore")                                                   # Ignore underflow error message
    plt.streamplot(XX,YY,Vx,Vy, linewidth=0.5, density=40, color='r',           # Plot streamlines
                    arrowstyle='-', start_points=XYsl)
    plt.clim(vmin=-0.51, vmax=2)
    plt.fill(xBoundary,yBoundary,'k')                                                         # Plot airfoil as black polygon
    plt.xlabel('X Units')                                                       # Set X-label
    plt.ylabel('Y Units')                                                       # Set Y-label
    plt.gca().set_aspect('equal')                                               # Set axes equal
    plt.xlim(xVals)                                                             # Set X-limits
    plt.ylim(yVals)                                                             # Set Y-limits
    plt.show()                                                                  # Display plot


    #tryk koefficient contour plot
    fig = plt.figure(3)                                                         # Create figure
    plt.cla()                                                                   # Get ready for plotting
    plt.contourf(XX,YY,pressureCoeffXY,500,cmap='jet')                                     # Plot contour
    #plt.clim(vmin=0, vmax=2)
    plt.fill(xBoundary,yBoundary,'k')                                                         # Plot airfoil as black polygon
    plt.xlabel('X Units')                                                       # Set X-label
    plt.ylabel('Y Units')                                                       # Set Y-label
    plt.gca().set_aspect('equal')                                               # Set axes equal
    plt.xlim(xVals)                                                             # Set X-limits
    plt.ylim(yVals)                                                             # Set Y-limits
    plt.show()

"""
   


   
   






