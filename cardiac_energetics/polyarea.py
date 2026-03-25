import numpy as np

def polyarea(x, y):
    """
    Calculate the area of a polygon using the Shoelace formula.
    
    Parameters:
    -----------
    x : array-like
        X coordinates of the polygon vertices
    y : array-like
        Y coordinates of the polygon vertices
        
    Returns:
    --------
    float
        Area of the polygon (always positive)
        
    Note:
    -----
    This is equivalent to MATLAB's polyarea function.
    Uses the Shoelace formula: A = 0.5 * |sum(x[i]*y[i+1] - x[i+1]*y[i])|
    """
    x = np.asarray(x)
    y = np.asarray(y)
    
    # Shoelace formula
    area = 0.5 * np.abs(np.dot(x, np.roll(y, -1)) - np.dot(y, np.roll(x, -1)))
    
    return area