'''
This script computes the convex hull efficiently for a uniformly random planar point set.

'''

import numpy as np
import matplotlib.pyplot as plt

def turning_angle(hull, pt):
    x_2prev = hull[:,-2]
    x_1prev = hull[:,-1]
    v1 = x_2prev - x_1prev
    v2 = pt - x_1prev
    # Technically we should add pi/2 to this, but we will be doing this to every column in remaining pts and taking a max, so it won't matter. 
    return np.arctan2(np.linalg.det(np.vstack((v1,v2))), np.dot(v1,v2)) # DOPE Lin Alg!

# Given a point set in 2D (2 x n np array), returns a 2 x m np array of pts in the hull. First row is x coordinates, and the second row is y coordinates. 
def jarvis_match(pts):
    if pts.shape[1] <= 2:
        return pts
    chull_col_inds = [np.argmin(pts[0,:])]
    #The leftmost point must be in the hull. Add auxilary pt to the hull to give starting orientation. 
    chull = pts[:, chull_col_inds[0]].reshape(2,1) + np.array([[0],[1]])
    chull = np.append(chull, pts[:, chull_col_inds].reshape(2,1), axis=1)
    
    while (chull.shape[1] < 3 or np.all(chull[:,1] != chull[:,-1], axis=0)) and chull.shape[1] <= pts.shape[1]:
        remaining_pts = np.delete(pts, chull_col_inds[-1], axis=1)
        turning = np.apply_along_axis(lambda a: turning_angle(chull,a), 0, remaining_pts)
        hull_new = remaining_pts[:,np.argmax(turning)]
        chull_col_inds.append(np.where(np.all(pts == hull_new[:,None], axis = 0))[0][0])
        chull = np.append(chull, hull_new[:,None],axis=1)

    return np.append(chull, chull[:,1][:,None], axis=1)[:,1:]
def graham_scan(pts):
    barycenter = np.mean(pts, axis=1)
    radii = np.linalg.norm(pts - np.tile(barycenter, (pts.shape[1],1)), axis = 1)
    theta = np.arctan2(pts[:,1]- barycenter[1], pts[:,0] - barycenter[0])
    inds = np.argsort(theta)    #Clean this sort up
    radii = radii[inds]
    theta = theta[inds]


if __name__ == '__main__':
    NUM_PTS=1500
    points = np.random.uniform(low=-5.0, high=5.0, size=(2,NUM_PTS))
    hull= jarvis_match(points)

    fig, ax = plt.subplots() 
    ax.set_aspect('equal')
    plt.scatter(points[0], points[1])
    plt.plot(hull[0], hull[1])

    plt.show()