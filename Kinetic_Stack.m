function out = Kinetic_Stack(DPdist)

%% Stackemup
% This function takes a degree of polymerization distribution as an input, in the form of an Nx2
% matrix:
% [DoP   # of chains]
% The output will be a histogram of the tie chain lengths

% as an initial example, we will use a distribution with average DP of 400
% (Mn = 67200) and a standard deviation of 50 (this is not the same as PDI)
% and assumes the chains follow a normal distribution, which is usually not
% the case

global DPs T kb
T = 298;
kb = 1.38E-23;
DPs = DPdist;
% pi_length = 1500;           % how many chains will we stack

Stack = [0 pick_pol()-1]; % stack is zero-indexed at the first monomer of the first chain to be picked
[n m] = size(Stack);

L= [];
count = 0;

while count<5
    if n>1
        Rates = get_rates(Stack,n) % Rates is [nx1]
        Cuts = make_bins(Rates); % cuts is also [nx1] ranging from 0 to 1
        process = choose_process(rand,Cuts);
        Stack = perform_process(Stack,process);
        L = [L;(length(Stack))];
        disp(L(end))
        if length(L)-5>0
            diff = L(end)-L(end-5);
            if abs(diff)<2
                count = count+1;
            else
                count = 0;
            end
        end
    else
        Stack=initiate(Stack);
    end
    [n, m] = size(Stack);
end

out = Stack;
stackplot(Stack)

end

function out = initiate(Stack)
%% Initiate
% Add a chain if there is only 1

DP = pick_pol();
OV = Stack(1,:);

hitstack = randi(OV(2)-OV(1)+1)+OV(1)-1; % where does new chain hit stack
hitchain = randi(DP);                    % which monomer of new chain hits stack

new_chain = [hitstack-hitchain+1, hitstack+DP-hitchain];
out = [Stack;new_chain];
end

function Cuts = make_bins(Rates)
%% Make Bins
% take nx1 [Rates] and turn it into nx1 [Cuts] ranging from 0 to 1
n = length(Rates);
ratesum = sum(Rates);
Cuts = zeros(n,1);
Cuts(1) = Rates(1)/ratesum;

for i = 2:n-1
    Cuts(i) = Rates(i)/ratesum+Cuts(i-1);
end
Cuts(n) = 1;
end

function process = choose_process(rando,Cuts)
%% Choose Process
% Take a random number and determine an integer that tells you which
% process to perform
above = find(Cuts > rando);
process = above(1);
end

%% PROCESSES
% The actual processes are contained here

function out = perform_process(Stack,p)
%% Perform Process
% take a selected process and the existing stack and DO THAT SHIT

if p==1
    out = add_front(Stack);
elseif p==2
    out = add_back(Stack);
elseif p==3
    out = det_front(Stack);
elseif p==4
    out = det_back(Stack);
end

end

function out = add_front(Stack)

new_chain = collide(Stack,length(Stack),length(Stack)-1);
out = [Stack; new_chain];

end

function out = add_back(Stack)

new_chain = collide(Stack,1,2);
out = [new_chain;Stack];

end

function out = collide(Stack,front,support)

%% collide
% takes a chain from DPs and collides a random
% monomer with the frontier chain (front, 1x2 vector) and any location
% where front overlaps with support (1x2 vector)

% the output is a 1x2 vector

DP = pick_pol();
[OL,OV] = find_overlap(Stack,front,support); % find overlap of previous two chains

hitstack = randi(OV(2)-OV(1)+1)+OV(1)-1; % where does new chain hit stack
hitchain = randi(DP);                    % which monomer of new chain hits stack

out = [hitstack-hitchain+1, hitstack+DP-hitchain];

end

function out = pick_pol()
global DPs
out = DPs(discr

etesample(DPs(:,3),1),1);
end

function out = det_front(Stack)
out = Stack(1:end-1,:);
end

function out = det_back(Stack)
out = Stack(2:end,:);
end

%% RATES
% Everything involved with getting rates is in this section
function out = get_rates(Stack,n)
% 1: add front
% 2: add back
% 3: detach front
% 4: detach back

processes = 4;
out = zeros(processes,1);

out(1) = add_front_rate(Stack,n);
out(2) = add_back_rate(Stack,n);
out(3) = det_front_rate(Stack,n);
out(4) = det_back_rate(Stack,n);

end

function out = add_front_rate(Stack,n)
k = n*exp(-n/100);
out = k;
end

function out = add_back_rate(Stack,n)
k = n*exp(-n/100);
out = k;
end

function out = det_front_rate(Stack,n)
global T kb
d0 = intrinsic_detach(n);
[overlap,OV] = find_overlap(Stack,length(Stack),length(Stack)-1);
Eov = overlap_E();
out = d0*exp(-overlap*Eov/(kb*T));
end

function out = det_back_rate(Stack,n)
global T kb
d0 = intrinsic_detach(n);
[overlap,OV] = find_overlap(Stack,1,2);
Eov = overlap_E();
out = d0*exp(-overlap*Eov/(kb*T));
end

function out = intrinsic_detach(n)
out = 8E11*exp(n/100);
end

function out = overlap_E()
out = 1.2E-21;
end


function [overlap_length,overlap_vector] = find_overlap(Stack,i,j)

overlap_vector = [max(Stack(i,1),Stack(j,1)), min(Stack(i,2),Stack(j,2))]; % find overlap of previous two chains
overlap_length = overlap_vector(2)-overlap_vector(1)+1;

end