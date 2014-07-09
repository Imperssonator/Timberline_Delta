function out = Mobility_Edge()

Ne = 10.^linspace(10,21,15)';

MOB = [Ne, zeros(15,4)];

for i = 2:5
    for j = 1:15
        MOB(j,i) = KMC(MOB(j,1),i);
    end
end

loglog(MOB(:,1),MOB(:,2:5),'or')

end

function out = KMC(Ne,i)

out = i;

end
