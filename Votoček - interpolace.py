import matplotlib.pyplot as plt
import numpy as np

while True: #ošetření případů, kdy uživatel zadá špatný název souboru nebo je soubor ve špatné formě
    try:
        nazev_souboru=input("Zadejte název souboru, ze kterého program načte data:")
        with open(nazev_souboru,"r") as body:
            x=[] #seznam x-ových souřadnic
            y=[] #seznam y-ových souřadnic
            for line in body:
                bod=line.split()
                x.append(float(bod[0]))
                y.append(float(bod[1]))
        break
    except FileNotFoundError:
        print("Soubor nenalezen")
    except ValueError:
        print("Špatný formát souboru")

"""seřazení bodů pole velikosti x-ových souřadnic"""
def quicksort(x,z=None): #x je seznam který chceme seřadit, z je podseznam indexů, na který je rekurzivně volána tato funkce
    if z==None:
        z=[i for i in range(len(x))]
    n=len(z)
    if n==1 or n==0:
        return z
    k=n//2
    a=[]
    b=[]
    for i in range(n):
        if x[z[i]]<x[z[k]]:
            a.append(z[i])
        if x[z[i]]>x[z[k]]:
            b.append(z[i])
    return quicksort(x,a) + [z[k]] + quicksort(x,b) #funkce vrací seznam indexů prvků původního seznamu podle velikosti

s=quicksort(x)
x=[x[i] for i in s]
y=[y[i] for i in s] #uspořádání vstupních dat

"""binární vyhledávání, bude se hodit pro výpočet chyby u některých metod"""
def binary_search(x,t,začátek,konec): #x je seznam, ve kterém hledáme, t je číslo, které hledáme
    k=(konec+začátek)//2
    if t>=x[k] and t<=x[k+1]:
        return k
    elif t<x[k]:
        return binary_search(x,t,začátek,k)
    else:
        return binary_search(x,t,k,konec) #funkce vrátí index prvku t, který je nejbližší nižší

w=[]

"""x,y jsou souřadnice zadaných bodů; w je seznam chyb jednotlivých metod, pokud uživatel chce nalézt metodu s nejmenší chybou; x1,y1 jsou
souřadnice zadaných bodů, podle kterých se chyba počítá"""
def linear(x,y,w=None,x1=None,y1=None):
    n=len(x)
    def linearni_interpolace(i,t): #i je index nejbližší menší hodnoty x ze seznamu, t je bod, ve kterém intepolujeme
        s=((y[i+1]-y[i])/(x[i+1]-x[i]))*(t-x[i])+y[i]
        return s
    def chyba(x1,y1): #výpočet chyby
        n=len(x)
        m=len(x1)
        chyba=0
        for i in range(m):
            k=linearni_interpolace(binary_search(x,x1[i],0,n-1),x1[i])
            chyba=chyba+(y1[i]-k)**2
        chyba=chyba/m
        return chyba
    if w!=None:
        w.append(chyba(x1,y1))
    if w==None: #vykreslení grafu je v podmínce, aby se při hledání nejmenší chyby nevykreslily grafy všech metod
        for i in range(n-1):
            a=np.linspace(x[i],x[i+1],50)
            b=linearni_interpolace(i,a)
            plt.plot(a,b)
        plt.scatter(x,y)
        plt.xlabel("x")
        plt.ylabel("y")
        plt.title("Lineární interpolace")
        plt.grid(True)
        plt.show()

def gaussova_eliminace(m): #Bude se hodit pro interpolaci polynomen a kubickými splajny
    n=len(m)
    for i in range (n):
        max=abs(m[i][i]) #Tady jsem se snažil o nějakou numerickou stabilitu
        k=i
        j=i
        while j<n:
            if abs(m[j][i])>max:
                max=abs(m[j][i])
                k=j
            j=j+1
        m[k],m[i]=m[i],m[k]
        l=m[i][i]
        for j in range (i,n+1):
            m[i][j]=m[i][j]/l
        for j in range(i+1,n):
            l=m[j][i]
            for k in range(i,n+1):
                m[j][k]=m[j][k]-l*m[i][k]
    p=[]
    for i in range(n-1,-1,-1):
        k=m[i][-1]/m[i][i]
        p.append(k)
        for j in range(n):
            m[j][-1]=m[j][-1]-k*m[j][i]
    return p #vrátí řešení v opačném pořadí

def lagrange(x,y):
    """Metoda interpolace polynomem pomocí n pomocných polynomů p(i), které jsou sestaveny tak, aby jejich hodnota v bodě i
    byla y[i] a v ostatních bodech x byla 0, součtem těchto n polynomů je interpolační polynom"""
    """Protože výsledná interpolační funkce by měla být až na drobné odchylky stejná u všech interpolací polynomem, zahrnul jsem výpočet
    chyby pouze ve funkci "lagrange_2" """
    a=np.linspace(x[0],x[-1],100)
    n=len(x)
    polynom=0
    for i in range(n):
        pomocny_polynom=y[i]
        for j in range(n):
            if i!=j:
                pomocny_polynom=pomocny_polynom*(a-x[j])/(x[i]-x[j])
        polynom=polynom+pomocny_polynom
    plt.plot(a,polynom)
    plt.scatter(x,y)
    plt.xlabel("x")
    plt.ylabel("y")
    plt.grid(True)
    plt.title("Lagrangeova interpolace")
    plt.show()
    
def lagrange_2(x,y,w=None,x1=None,y1=None):
    n=len(x)
    """Hledám polynom p stupně n-1, který splňuje p(x[i]))=y[i] pro každé i mezi 0 a n, to vede na gaussovu eliminaci"""
    m=[]
    for i in range(n):
        a=[1]
        for j in range(n-1):
            a.append(a[-1]*x[i])
        a.append(y[i])
        m.append(a)
    p=gaussova_eliminace(m) #koeficienty polynomu
    a=np.linspace(x[0],x[-1],100)
    def lagrangeuv_polynom(p,t): #p je seznam koeficientů polynomu, t je bod, ve kterém interpolujeme
        b=p[0]
        m=len(p)
        for i in range(1,m):
            b=b*t+p[i] #Hornerovo schéma
        return b
    def chyba(x1,y1):
        m=len(x1)
        chyba=0
        for i in range(m):
            k=lagrangeuv_polynom(p,x1[i])
            chyba=chyba+(y1[i]-k)**2
        chyba=chyba/m
        return chyba
    if w!=None:
        w.append(chyba(x1,y1))
    b=lagrangeuv_polynom(p,a)
    if w==None:
        plt.plot(a,b)
        plt.scatter(x,y)
        plt.xlabel("x")
        plt.ylabel("y")
        plt.grid(True)
        plt.title("Lagrangeova interpolace")
        plt.show()

def newton(x,y):
    """třetí metoda interpolace polynomem, koeficienty získám pomocí poměrných diferencí"""
    n=len(x)
    from collections import deque
    f=deque([]) #využívám frontu
    for i in range(n):
        f.append(y[i])
    s=[y[0]] #koeficienty polynomu
    for i in range(1,n):
        for j in range(i,n):
            a=f.popleft()
            b=f[0]
            f.append((b-a)/(x[j]-x[j-i]))
        f.popleft()
        s.append(f[0])
    def newtonuv_polynom(s,t): #s je seznam koeficientů, t je bod, ve kterém interpolujeme
        n=len(s)
        c=0
        for i in range(n):
            b=s[i]
            for j in range(i):
                b=b*(t-x[j]) #Hornerovo schéma
            c=c+b
        return c
    a=np.linspace(x[0],x[-1],100)
    b=newtonuv_polynom(s,a)
    plt.plot(a,b)
    plt.scatter(x,y)
    plt.xlabel("x")
    plt.ylabel("y")
    plt.grid(True)
    plt.title("Newtonova interpolace")
    plt.show()

def newton_2(x,y,w=None,x1=None,y1=None):
    """Newtonova metoda má tu výhodu, že je oproti jiným metodám výpočetně méně náročné přidat další body. Pokusil jsem se přepsat původní
    funkci "newton", abych toho mohl využít. Vypozoroval jsem, že výsledek se celkem znatelně liší, tak jsem v programu nechal
    obě verze."""
    n=len(x)
    m=[[0]*n for i in range(n)] #tabulka poměrných diferencí
    for i in range(n):
        m[i][0]=y[i]
    for i in range(1,n):
        for j in range(i,n):
            m[j][i]=(m[j][i-1]-m[j-1][i-1])/(x[j]-x[j-i])
    def newtonuv_polynom(m,t): #m je tabulka poměrných diferencí, koeficienty polynomu jsou na diagonále, t je interpolační bod
        c=m[-1][-1]
        for i in range(n-1,-1,-1):
            c=c*(t-x[i])
            c=c+m[i][i] #Hornerovo schéma
        return c
    def chyba(x1,y1):
        l=len(x1)
        chyba=0
        for i in range(l):
            k=newtonuv_polynom(m,x1[i])
            chyba=chyba+(y1[i]-k)**2
        chyba=chyba/l
        return chyba
    if w!=None:
        w.append(chyba(x1,y1))
    if w==None:
        a=np.linspace(x[0],x[-1],100)
        b=newtonuv_polynom(m,a)
        plt.plot(a,b)
        plt.scatter(x,y)
        plt.xlabel("x")
        plt.ylabel("y")
        plt.grid(True)
        plt.title("Newtonova interpolace")
        plt.show()
    return m

"""Možnost přidat další bod do Newtonovy interpolace. m je původní tabulka; x,y seznamy souřadnic všech bodů; x2,y2 jsou souřadnice přidávaného
bodu. Abych mohl přidat více bodů najednou, přidal jsem ty podmínky "if x2!=None" a "if x2==None". Podle toho program pozná, že už všechny body 
přidal a vykreslí graf"""
def dalsi_bod(m,x,y,x2,y2):
    n=len(m)
    x.append(x2)
    y.append(y2)
    """rozšíření tabulky"""
    m.append([0]*(n+1))
    m[-1][0]=y2
    if x2!=None: #Pokud x==None, znamená to, že už všechny body byly přidány
        for i in range(1,n+1):
            m[-1][i]=(m[-1][i-1]-m[-2][i-1])/(x2-x[-i-1])
    def newtonuv_polynom(m,t):
        c=m[-1][-1]
        for i in range(n-1,-1,-1):
            c=c*(t-x[i])
            c=c+m[i][i]
        return c
    if x2==None:
        a=np.linspace(min(x[0:-1]),max(x[0:-1]),100)
        b=newtonuv_polynom(m,a)
        plt.plot(a,b)
        plt.scatter(x,y)
        plt.xlabel("x")
        plt.ylabel("y")
        plt.grid(True)
        plt.title("Newtonova interpolace")
        plt.show()
    return m

def quadratic_spline(x,y,w=None,x1=None,y1=None): #analogie kubických splajnů
    n=len(x)
    m=[[0]*((n-1)*3+1) for i in range((n-1)*3)]
    for i in range(n-2):
        h=x[i+1]-x[i]
        m[i*3][i*3+2]=1
        m[i*3][-1]=y[i]
        m[i*3+1][i*3]=h**2
        m[i*3+1][i*3+1]=h
        m[i*3+1][-1]=y[i+1]-y[i]
        m[i*3+2][i*3]=2*h
        m[i*3+2][i*3+1]=1
        m[i*3+2][i*3+4]=-1
    m[-3][-2]=1
    m[-3][-1]=y[-2]
    h=x[-1]-x[-2]
    m[-2][-4]=h**2
    m[-2][-3]=h
    m[-2][-1]=y[-1]-y[-2]
    m[-1][-4]=-1
    m[-1][-7]=1
    p=gaussova_eliminace(m)
    def interpolacni_funkce(t,p,x,i):
        s=p[-1-i*3]*(t-x[i])**2+p[-2-i*3]*(t-x[i])+p[-(i+1)*3]
        return s
    a=np.linspace(x[0],x[-1],50)
    def chyba(x1,y1):
        m=len(x1)
        chyba=0
        for i in range(m):
            k=interpolacni_funkce(x1[i],p,x,binary_search(x,x1[i],0,n-1))
            chyba=chyba+(y1[i]-k)**2
        chyba=chyba/m
        return chyba
    if w!=None:
        w.append(chyba(x1,y1))
    if w==None:
        for i in range(n-1):
            a=np.linspace(x[i],x[i+1],50)
            b=interpolacni_funkce(a,p,x,i)
            plt.plot(a,b)
        plt.scatter(x,y)
        plt.xlabel("x")
        plt.ylabel("y")
        plt.grid(True)
        plt.title("Kvadratické splajny")
        plt.show()
    

def cubic_spline(x,y,w=None,x1=None,y2=None):
    """Mezi každými dvěma body proložím kubický polynom p[i]= a_i*(x-x_i)**3+b_i*(x-x_i)**2+c_i*(x-x_i)+d_i, tak, aby byly splněny podmínky:
    1) p(x[i])=y[i]
    2) p(x[i+1])=y[i+1]
    3) derivace v bodě [i+1] funkce p[i] se rovná derivaci funkce p[i+1] ve stejném bodě
    4) rovnají se i druhé derivace
    Vznikne soustava rovnic, řešením jsou koeficienty polynomů"""
    n=len(x)
    m=[[0]*((n-1)*4+1) for i in range((n-1)*4)]
    for i in range(n-2):
        h=x[i+1]-x[i]
        """podmínka 1"""
        m[i*4+1][i*4+3]=1
        m[i*4+1][-1]=y[i]
        """podmínka 2"""
        m[i*4+2][i*4]=h**3
        m[i*4+2][i*4+1]=h**2
        m[i*4+2][i*4+2]=h
        m[i*4+2][-1]=y[i+1]-y[i]
        """podmínka 3"""
        m[i*4+3][i*4]=3*h**2
        m[i*4+3][i*4+1]=2*h
        m[i*4+3][i*4+2]=1
        m[i*4+3][i*4+6]=-1
        """podmínka 4"""
        m[(i+1)*4][i*4]=6*h
        m[(i+1)*4][i*4+1]=2
        m[(i+1)*4][i*4+5]=-2
    """podmínky 1 a 2 pro poslední bod"""
    h=x[-1]-x[-2]
    m[-3][-2]=1
    m[-3][-1]=y[-2]
    m[-2][-5]=h**3
    m[-2][-4]=h**2
    m[-2][-3]=h
    m[-2][-1]=y[-1]-y[-2]
    """Další 2 podmínky navíc, aby matice byla regulární. Zajistí to, že první dva polynomy jsou stejné, stejně tak poslední dva"""
    m[0][0]=1
    m[0][4]=-1
    m[-1][-5]=1
    m[-1][-9]=-1
    p=gaussova_eliminace(m)
    def interpolacni_funkce(t,p,x,i):
        s=p[(n-2-i)*4+3]*(t-x[i])**3+p[(n-2-i)*4+2]*(t-x[i])**2+p[(n-2-i)*4+1]*(t-x[i])+p[(n-2-i)*4]
        return s
    def chyba(x1,y1):
        m=len(x1)
        chyba=0
        for i in range(m):
            k=interpolacni_funkce(x1[i],p,x,binary_search(x,x1[i],0,n-1))
            chyba=chyba+(y1[i]-k)**2
        chyba=chyba/m
        return chyba
    if w!=None:
        w.append(chyba(x1,y1))
    if w==None:
        for i in range(n-1):
            a=np.linspace(x[i],x[i+1])
            b=interpolacni_funkce(a,p,x,i)
            plt.plot(a,b)
        plt.scatter(x,y)
        plt.xlabel("x")
        plt.ylabel("y")
        plt.grid(True)
        plt.title("Kubické spliny")
        plt.show()

"""Primitivní metoda interpolace, každému bodu přiřadí funkční hodnotu nejbližšího zadaného bodu"""
def nearest_neighbour(x,y,w=None,x1=None,y2=None):
    n=len(x)
    s=[] #seznam "hraničních" hodnot
    for i in range(n-1):
        s.append((x[i]+x[i+1])/2)
    if w==None:
        a=np.linspace(x[0],s[0],10)
        b=[y[0]]*10
        plt.plot(a,b)
        for i in range(n-2):
            a=np.linspace(s[i],s[i+1],10)
            b=[y[i+1]]*10
            plt.plot(a,b)
        a=np.linspace(s[-1],x[-1],10)
        b=[y[-1]]*10
        plt.plot(a,b)
    def chyba(x1,y1):
        m=len(x1)
        chyba=0
        for i in range(m):
            if x1[i]<s[0]:
                k=x[0]
            elif x1[i]>s[-1]:
                k=x[-1]
            else:
                k=binary_search(s,x1[i],0,m-1)
            chyba=chyba+(y1[i]-y[k+1])**2
        chyba=chyba/m
        return chyba
    if w!=None:
        w.append(chyba(x1,y1))
    if w==None:
        plt.scatter(x,y)
        plt.xlabel("x")
        plt.ylabel("y")
        plt.grid(True)
        plt.title("Metoda nejbližšího souseda")
        plt.show()

"Metoda interpolace polynomem, k výpočtu se interpolované hodnoty se používají váhy závislé na vzdálenosti tohoto bodu od interpolačních bodů"
def barycentricka_interpolace(x,y):
    n=len(x)
    s=[] #seznam vah
    for i in range(n): #výpočet vah
        k=1
        for j in range(n):
            if i!=j:
                k=k/(x[i]-x[j])
        s.append(k)
    def interpolacni_funkce(x,y,s,t): #x,y jsou souřadnice zadaných bodů, s je seznam vah, t je interpolační bod
        l=0
        for i in range(n):
            try:
                l=l+s[i]*y[i]/(t-x[i])
            except ZeroDivisionError: #nulou dělím právě tehdy když t==x[i], v takovém případě chci, aby funkce vrátila y[i]
                l=y[i]
        m=0
        for i in range(n):
            try:
                m=m+s[i]/(t-x[i])
            except ZeroDivisionError:
                m=1
        return l/m
    a=np.linspace(x[0],x[-1],100)
    b=interpolacni_funkce(x,y,s,a)
    plt.plot(a,b)
    plt.scatter(x,y)
    plt.xlabel("x")
    plt.ylabel("y")
    plt.grid(True)
    plt.title("Barycentrická interpolace")
    plt.show()

"""typ interpolace, který si uživatel přeje, pokud si vybere "nejmenší chyba", program provede interpolaci s nejmenší naměřenou chybou na
základě dat, která mu uživatel poskytne"""
typ_interpolace=input("Zvolte typ interpolace (lineární, polynomiální, kvadratické splajny, kubické splajny, metoda nejbližšího souseda, nejmenší chyba):")
s=["lineární", "polynomiální", "kvadratické splajny", "kubické splajny", "metoda nejbližšího souseda", "barycentrická", "nejmenší chyba"]
def preklep(typ_interpolace,s): #Tady jsem ošetřil případ, kdy se uživatel přepíše
    if typ_interpolace in s:
        return typ_interpolace
    else:
        typ_interpolace=typ_interpolace.lower()
        import difflib
        typ_interpolace=difflib.get_close_matches(typ_interpolace,s)
        if typ_interpolace:
            return typ_interpolace[0]
        else:
            typ_interpolace=input("Typ interpolace nebyl rozpoznán, vyberte, prosím, znovu z nabízených možností: (lineární, polynomiální, kvadratické splajny, kubické splajny, metoda nejbližšího souseda, nejmenší chyba)")
            return preklep(typ_interpolace,s)
def preklep_2(metoda,s):
        if metoda in s:
            return metoda
        else:
            metoda=metoda.lower()
            import difflib
            metoda=difflib.get_close_matches(metoda,s)
            if metoda:
                return metoda[0]
            else:
                metoda=input("Metoda polynomiální interpolace nebyla rozpoznána, vyberte, prosím, znovu z nabízených možností: (lagrange, soustava rovnic, newton, barycentrická)")
                return preklep_2(metoda,s)

typ_interpolace=preklep(typ_interpolace,s)

if typ_interpolace=="lineární":
    linear(x,y)
elif typ_interpolace=="polynomiální":
    metoda=input("Zvolte metodu polynomiální interpolace (lagrange, soustava rovnic, newton, newton 2, barycentrická):")
    s=["lagrange", "soustava rovnic", "newton","newton 2", "barycentrická"]
    if preklep_2(metoda,s)=="lagrange":
        lagrange(x,y)
    elif preklep_2(metoda,s)=="soustava rovnic":
        lagrange_2(x,y)
    elif preklep_2(metoda,s)=="newton":
        newton(x,y)
    elif preklep_2(metoda,s)=="newton 2":
        m=newton_2(x,y)
        dalsibody=input("Přejete si přidat další body? (ano, ne)")
        while dalsibody!="ne" and dalsibody!="ano":
            dalsibody=input("Odpověď nerozpoznána. Přejete si přidat další body?(ano,ne):")
        if dalsibody == "ne":
            pass
        else:
            while True:
                try:
                    nazev_souboru=input("Zadejte název souboru:")
                    with open(nazev_souboru,"r") as db:
                        x2=[]
                        y2=[]
                        for line in db:                       
                            bod=line.split()
                            x2.append(float(bod[0]))
                            y2.append(float(bod[1]))
                        x2.append(None)
                        y2.append(None)
                        n=len(x2)
                        for i in range(n):
                            m=dalsi_bod(m,x,y,x2[i],y2[i])
                    break
                except FileNotFoundError:
                    print("Soubor nenalezen")
                except ValueError:
                    print("Špatný formát souboru")
    else:
        barycentricka_interpolace(x,y)
elif typ_interpolace=="kvadratické splajny":
    quadratic_spline(x,y)
elif typ_interpolace=="kubické splajny":
    cubic_spline(x,y)
elif typ_interpolace=="metoda nejbližšího souseda":
    nearest_neighbour(x,y)
else:
    while True:
        try:
            nazev_souboru=input("Zadejte název souboru:")
            with open(nazev_souboru,"r") as kon:
                x1=[]
                y1=[]
                for line in kon:
                    bod=line.split()
                    x1.append(float(bod[0]))
                    y1.append(float(bod[1]))
            break
        except FileNotFoundError:
            print("Soubor nenalezen")
        except ValueError:
            print("Špatný formát souboru")
    linear(x,y,w,x1,y1)
    lagrange_2(x,y,w,x1,y1)
    newton_2(x,y,w,x1,y1)
    quadratic_spline(x,y,w,x1,y1)
    cubic_spline(x,y,w,x1,y1)
    nearest_neighbour(x,y,w,x1,y1)
    k=0
    for i in range(1,4):
        if w[i]<w[k]:
            k=i
    if k==0:
        print("Interpolace s nejmenší naměřenou chybou je lineární interpolace")
        linear(x,y)
    elif k==1:
        metoda=input("Interpolace s nejmenší naměřenou chybou je polynomiální interpolace, vyberte metodu (lagrange, soustava rovnic, newton, barycentrická):")
        metoda=preklep_2(metoda,s)
        if metoda=="lagrange":
            lagrange(x,y)
        elif metoda=="soustava rovnic":
            lagrange_2(x,y)
        elif metoda=="newton":
            newton(x,y)
        else:
            barycentricka_interpolace(x,y)
    elif k==2:
        print("Interpolace s nejmenší naměřenou chybou je Newtonova interpolace")
        newton_2(x,y)
    elif k==3:
        print("Interpolace s nejmenší naměřenou chybou je interpolace kvadratickými splajny")
        quadratic_spline(x,y)
    elif k==4:
        print("Interpolace s nejmenší naměřenou chybou je interpolace kubickými splajny")
        cubic_spline(x,y)
    else:
        print("Interpolace s nejmenší naměřenou chybou je barycentrická interpolace")
        barycentricka_interpolace(x,y)
