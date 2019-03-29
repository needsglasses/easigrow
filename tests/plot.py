from pylab import *
import numpy as np

ion()

x = np.genfromtxt("test-comparison",skip_header=1,names=True)

figure()
title('Easigrow comparison using --dadn white:barter14-aa7050t7451')
plot( x['test'][x['coupon']==0], x['easigrow'][x['coupon']==0],'sr',label='Short crack')
plot( x['test'][x['coupon']==1], x['easigrow'][x['coupon']==1],'vb',label='Middle tension')
plot( x['test'][x['coupon']==2], x['easigrow'][x['coupon']==2],'og',label='Compact')
plot([0,600], [0,600],'--')
xlim(0,600)
ylim(0,600)
xlabel('Test (blocks)')
ylabel('Easigrow (blocks)')
legend()

savefig('test-comparison.svg')
