

def Dint(arr,h):
    #复合辛普森积分
    n=len(arr)
    sum=0
    if (n%2==0):
        sum+=(arr[-1]+arr[-2])*h/2
        n-=1
    tsum=arr[0]+arr[n-1]
    for i in range(1,n-1):
        if(i%2):
            tsum+=4*arr[i]
        else:
            tsum+=2*arr[i]
    sum+=tsum*h/3
    return sum