import sys

a,da,b,db,at,dat,bt,dbt,ap,dap,bp,dbp = [float(var) for var in sys.argv[1::]]

alpha = a*ap/at
dalpha = ( (da*ap/at)**2 + (a*dap/at)**2 + (a*ap*dat/at**2)**2 )**.5

beta = bp + ap*b - a*ap*bt/at
#dbeta = ( (dbp)**2 + (dap*b)**2 + (ap*db)**2 + (da*ap*bt/at)**2 + (a*dap*bt/at)**2 + (a*ap*dbt/at)**2 + (a*ap*bt*dat/at**2)**2 )**.5
dbeta = ( (dbp)**2 + (dap*(b-a*bt/at))**2 + (db*ap)**2 + (da*ap*bt/at)**2 + (dbt*a*ap/at)**2 + (dat*a*ap*bt/at**2)**2 )**.5

print('a =% .2f (%.2f), b =% .1f (%.1f)'%(a,da,b,db))
print('at=% .2f (%.2f), bt=% .1f (%.1f)'%(at,dat,bt,dbt))
print('ap=% .2f (%.2f), bp=% .1f (%.1f)'%(ap,dap,bp,dbp))
print('Result:')
print('a =% .2f (%.2f), b =% .1f (%.1f)'%(alpha,dalpha,beta,dbeta))
