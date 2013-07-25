module flow
using PyCall
@pyimport scipy.spatial as sp


#Need to get a 6x6 viscosity matrix C
#routine to build sparse matrix by COLUMN
#to solve 
#[u_x]' [
#[u_y]  [
#[u_z]  [
#[p]    [

