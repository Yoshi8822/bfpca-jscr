library(reticulate)
library(introdataviz)
library(fda)

# Import the necessary Python modules
pickle <- import("pickle")
builtins <- import("builtins")
numpy <- import("numpy")

# Load the pickle file in R
py_array <- builtins$open("This is where I upload a pkl file")
r_array <- pickle$load(py_array)

# Convert the Python array to an R array
r_array <- numpy$array(r_array)

# Check the shape of the R array
dim(r_array)
print(r_array[,,1])


##################  Use the roughness penalty to determine basis coefficients ###############

# Create a B-spline basis for the smoothing
bspl_40 <- create.bspline.basis(c(0, 100), 40, 4)
fdParobj_pols <- fdPar(fdobj = bspl_40, Lfdobj = int2Lfd(m = 2), lambda = 1)

# Apply smoothing to the multidimensional r_array
fdatime = seq(0, 100, len= 101)
pls_fdSmooth <- smooth.basis(fdatime, r_array, fdParobj = fdParobj_pols)$fd

# Plot the original data points for the first sample
plot(r_array[, 1, 1], r_array[, 1, 2], pch = 20, 
     xlab = "Horizontal Displacement", ylab = "Vertical Displacement")

# Overlay the smoothed curves with different lambdas
log10_lambda_range <- seq(-10, 6, by = 2)
lambda_range <- 10^log10_lambda_range
n_lambda <- length(lambda_range)

for(lam in seq_len(n_lambda)) {
  fdParobj_lam <- fdPar(fdobj = bspl_40, 
                        Lfdobj = int2Lfd(m = 2), 
                        lambda = lambda_range[lam])
  
  # Smooth the multidimensional array with the current lambda
  pls_fdSmooth_lam <- smooth.basis(fdatime, r_array, fdParobj = fdParobj_lam)$fd
  # Extract smoothed data for the first sample
  y_smooth_lam <- eval.fd(fdatime, pls_fdSmooth_lam[,1])
  z_smooth_lam <- eval.fd(fdatime, pls_fdSmooth_lam[,2])
  
  # Overlay smoothed trajectory for the first sample
  lines(y_smooth_lam, z_smooth_lam, col = lam)
}

legend("topright", legend = paste0("lambda = ", lambda_range), col = seq_len(n_lambda), lty = 1)

# Assuming r_array is structured as r_array[time_points, samples, dimensions]
gcv_vec <- vector(mode = "numeric", length = n_lambda)

for (lam in seq_len(n_lambda)) {
  # Create fdPar with chosen lambda
  fdParobj_lam <- fdPar(fdobj = bspl_40, Lfdobj = int2Lfd(m = 2), lambda = lambda_range[lam])
  
  # Apply smoothing to the multidimensional array
  pls_fdSmooth_lam <- smooth.basis(argvals = fdatime, y = r_array, fdParobj = fdParobj_lam)
  
  # Calculate Generalised Cross Validation (GCV)
  gcv_dim1 <- mean(pls_fdSmooth_lam$gcv[, 1])
  gcv_dim2 <- mean(pls_fdSmooth_lam$gcv[, 2])
  gcv_vec[lam] <- mean(c(gcv_dim1, gcv_dim2))  # Average GCV over dimensions
}

# Plot GCV values across log10(lambda) range
plot(log10_lambda_range, 
     y = gcv_vec, 
     type = "b", 
     xlab = expression(log[10](lambda)), 
     ylab = "GCV")

# Identify the best lambda based on minimum GCV
best_lambda_index <- which.min(gcv_vec)
abline(v = log10_lambda_range[best_lambda_index])

# Set up final smoothed object
fdParobj_final <- fdPar(fdobj = bspl_40, 
                        Lfdobj = int2Lfd(m = 2),
                        lambda = lambda_range[best_lambda_index]) # create fdPar with chosen 
pls_fdSmooth_final <- smooth.basis(argvals = fdatime, y = r_array, fdParobj = fdParobj_final)
final_fd <- pls_fdSmooth_final$fd



################## Carry out the PCA of the bivariate functional data object 'final_fd' using three harmonics ##################
nharm = 3
fdapcaList = pca.fd(final_fd, nharm)

# Perform a VARIMAX rotation of these eigenfunctions to make it easier to interpret
fdarotpcaList = varmx.pca.fd(fdapcaList)

# The following code can be used to determine the number of harmonics used
fdaeig = fdapcaList$values
neig = 12
x = matrix(1,neig-nharm,2)
x[,2] = (nharm+1):neig
y = as.matrix(log10(fdaeig[(nharm+1):neig]))
c = lsfit(x,y,int=FALSE)$coef
op <- par(mfrow=c(1,1),cex=1.2)
plot(1:neig, log10(fdaeig[1:neig]), "b",
     xlab="Eigenvalue Number",
     ylab="Log10 Eigenvalue")
lines(1:neig, c[1]+ c[2]*(1:neig), lty=2)
par(op)
