function [out, UsedChains] = Kinetic_Stack(DPdist)

clc
%% Stackemup
% This function takes a degree of polymerization distribution as an input, in the form of an Nx2
% matrix:
% [DoP   # of chains]
% The output will be a histogram of the tie chain lengths

% as an initial example, we will use a distribution with average DP of 400
% (Mn = 67200) and a standard deviation of 50 (this is not the same as PDI)
% and assumes the chains follow a normal distribution, which is usually not
% the case

%% Hard Coded Parameters
global T kb DPdist1 Stack
DPdist1 = DPdist;
T = 298;
kb = 1.38E-23;
Iterations = 100000;
% pi_length = 300;           % how many chains will we stack

Stack = [0 DPSample()-1]; % stack is zero-indexed at the first monomer of the first chain to be picked

[m, n, p] = size(Stack);

L = ones(Iterations,1);
iter = 0;

while iter<Iterations 
    iter = iter+1;
    ProcessList = find_process();  % Number of possible processes is listed
    % 1 add chain to pi stack
    % 2 remove chain from pi stack
    % 3 new chain in alkyl stacking direction
    % 4 remove a single chain in the alkyl stack
    % 5 add a chain that forms both pi and alkyle stack
    % 6 remove a chain in both pi and alkyl stack
    Rates = get_rates(ProcessList); % Rates is [nx1]
    Cuts = make_bins(Rates); % cuts is also [nx1] ranging from 0 to 1
    process = choose_process(rand,Cuts);
    Stack = perform_process(process, ProcessList);
    disp(length(Stack))
    L(iter) = length(Stack);
    [m, n p] = size(Stack);
    
%     figure(1)
%     stackplot(Stack)
    
%     if mod(iter,10000) == 0
%         figure(1)
%         hold on
%         stackplot(Stack,iter)
%         drawnow
        
%         figure
%         hist(DPdist1,20);           % Histogram of remaining free chains
        
    end
end

% StartDist = tabulate(DPdist);
% First = StartDist(1,1); Last = StartDist(end,1);
% EndDist = tabulate(DPdist1);
% InnerFirst = EndDist(1,1); InnerLast = EndDist(end,1);
% Prefix = (First:InnerFirst-1)';
% Suffix = (InnerLast+1:Last)';
% Prefix = [Prefix, zeros(size(Prefix,1),2)];
% Suffix = [Suffix, zeros(size(Suffix,1),2)];
% EndDist = [Prefix; EndDist; Suffix];                % Add zeros for the lengths that were totally depleted in the simulation
% 
% UsedChains = [StartDist(:,1), StartDist(:,2:3)-EndDist(:,2:3)];
% figure
% bar(UsedChains(:,1),UsedChains(:,2))
% 
% figure(2)
% plot((1:length(L))',L,'-b')
% out = Stack;



end

function out = find_process()
global Stack occupied_top occupied_bottom growth_sites removables
[m, ~, p] = size(Stack);
out = [1 0; 2 0; 3 0; 4 0; 5 0; 6 0];
out(1,2) = 2;   % Two pi chains can be added in either direction

if m>1
    out(2,2) = 2;   % Two possible ways of removing in either direction
end

top_layer = Stack(:,:,1);
occupied_top = find_sites(top_layer);
bottom_layer = Stack(:,:,end);
occupied_bottom = find_sites(bottom_layer);
out(3,2) = size(occupied_top) + size(occupied_bottom);  % Number of places where new alkyl chain can begin

if size(occupied_top) == 1 && size(occupied_bottom) == 1
    out(4,2) = 2;
elseif size(occupied_top) == 1 || size(occupied_bottom) == 1
    out(4,2) = 1;
end
% Number of ways in which single alkyl stack chain can be removed

FirstLayer = locate_stack();
growth_sites = [];
removables = [];
for i = 1:FirstLayer-1
    occupied = find_sites(Stack(:,:,i));
    if occupied(1) ~= 1 && Stack(occupied(1)-1,:,i+1) ~= [0 0]
        out(5,2) = out(5,2) + 1;
        growth_sites = [growth_sites; [occupied(1) i 1]];
    end
    if occupied(end) ~= m && Stack(occupied(end)+1,:,i+1) ~= [0 0]
        out(5,2) = out(5,2) + 1;
        growth_sites = [growth_sites; [occupied(end) i 2]];
    end
    if i == 1 && size(occupied)>1
        out(6,2) = out(6,2)+2;
        removables = [removables; [occupied(1),i]; [occupied(end),i]];
    end
    if i>1 && Stack(occupied(end),:,i-1) == [0 0]
        out(6,2) = out(6,2) + 1;
        removables = [removables; [occupied(end),i]];
    end
    if i>1 && Stack(occupied(1),:,i-1) == [0 0]
        out(6,2) = out(6,2) + 1;
        removables = [removables; [occupied(1),i]];
    end
end
for i = p:-1:FirstLayer+1
    occupied = find_sites(Stack(:,:,i));
    if occupied(1) ~= 1 && Stack(occupied(1)-1,:,p-1) ~= [0 0]
        out(5,2) = out(5,2) + 1;
        growth_sites = [growth_sites; [occupied(1) i 2]];
    end
    if occupied(end) ~= m && Stack(occupied(end)+1,:,p-1) ~= [0 0]
        out(5,2) = out(5,2) + 1;
        growth_sites = [growth_sites; [occupied(end) i 2]];
    end
    if i == p && size(occupied) > 1
        out(6,2) = out(6,2)+2;
        removables = [removables; [occupied(1) i]; [occupied(end),i]];
    end
    if i<p && Stack(occupied(end),:,i+1) == [0 0]
        out(6,2) = out(6,2) + 1;
        removables = [removables; [occupied(end),i]];
    end
    if i<p && Stack(occupied(1),:,i+1) == [0 0]
        out(6,2) = out(6,2) + 1;
        removables = [removables; [occupied(1),i]];
    end
        
end
% Number of ways in which pi-alkyl chain can be added or removed
end

function out = locate_stack()
global Stack
[m, ~, ~] = size(Stack);
for j = 1:size(Stack,3)
    occupied = find_sites(Stack(:,:,j));
    if size(occupied) == m
        out = j;
        break
    end
end
end

function out = find_sites(single_layer)
find_zero = single_layer==0;    % Assigns zeros in the place of non zero elements and 1s in the place of zeros
find_nonzero = ~all(find_zero,2);   % Vector with 1s if the corresponding row in find_zero has at least one zero and 0 otherwise
out = find(find_nonzero);   % Indices of rows which have at least one non zero element ie. occupied sites
end

function out = initiate(Stack)
%% Initiate
% Add a chain if there is only 1

DP = DPSample();
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

function out = perform_process(p,List)
%% Perform Process
% take a selected process and the existing stack and DO THAT SHIT
global occupied_top occupied_bottom growth_sites
q = randi(List(p,2));
if p == 1
    if q == 1
        out = add_front();
    elseif q == 2
        out = add_back();
    end
elseif p == 2
    if q == 1
        out = det_front();
    elseif q == 2
        out = det_back();
    end
elseif p == 3
    if q<size(occupied_top)
        out = add_top(q);
    else
        out = add_bottom(q-size(occupied_top));
    end
elseif p == 4
    if List(p,2) == 1
        if size(occupied_top) == 1
            out = det_top();
        elseif size(occupied_bottom) == 1
            out = det_bottom();
        end
    else
        if q == 1
            out = det_top();
        elseif q == 2
            out = det_bottom();
        end
    end
elseif p == 5
    out = add_pialkyl(growth_sites(q,:));
elseif p == 6
    out = det_pialkyl(q);
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

function out = collide(layer,front,support)
%% collide
% takes a chain from DPs and collides a random
% monomer with the frontier chain (front, 1x2 vector) and any location
% where front overlaps with support (1x2 vector)

% the output is a 1x2 vector

DP = DPSample();
[OL,OV] = find_overlap(front,support,layer); % find overlap of previous two chains

hitstack = randi(OV(2)-OV(1)+1)+OV(1)-1; % where does new chain hit stack
hitchain = randi(DP);                    % which monomer of new chain hits stack

out = [hitstack-hitchain+1, hitstack+DP-hitchain];

end

function out = collide_alkyl(layer,chain)
global Stack
[m, n, p] = size(Stack);
out = zeros(m,n);
DP = DPSample();
hitstack = randi(Stack(chain,2,layer) - Stack(chain,1,layer) + 1);
hitchain = randi(DP);

new_chain = [hitstack-hitchain+1, hitstack+DP-hitchain];
out(chain,:) = new_chain;
end

function out = collide_pialkyl(front_chain)
global Stack 
DP = DPSample();
hitstack = randi(Stack(front_chain(1),2,front_chain(2))-Stack(front_chain(1),1,front_chain(2))+1);
hitchain = randi(DP);

out = [hitstack-hitchain+1, hitstack+DP-hitchain];

end

% function out = pick_pol()
% global DPs
% out = DPs(discretesample(DPs(:,3),1),1);
% end

function out = det_front()
global Stack
redefine_det(Stack(end,:));
out = Stack(1:end-1,:);

end

function out = det_back()
global Stack
redefine_det(Stack(1,:));
out = Stack(2:end,:);

end

function out = add_top(q)
global Stack
Stack = cat(3, Stack, collide_alkyl(q,1));
out = Stack;
end

function out = add_bottom(q)
global Stack
Stack = cat(3, collide_alkyl(q,size(Stack,3)));
out  = Stack;
end

function out = det_top()
global Stack
Stack(:,:,1) = [];
out = Stack;
end

function out = det_bottom()
global Stack
Stack(:,:,end) = [];
out = Stack;
end

function out = add_pialkyl(front_chain)
global Stack
new_chain = collide_pialkyl(front_chain);
if front_chain(3) == 1
    Stack(front_chain(1)-1,:,front_chain(2)) = new_chain;
else
    Stack(front_chain(1)+1,:,front_chain(2)) = new_chain;
end
out = Stack;
end

function out = det_pialkyl(q)
global Stack removables 
det = removables(q,:);
Stack(det(1),:,det(2)) = [0 0];
out =  Stack;
end

%% RATES
% Everything involved with getting rates is in this section
function out = get_rates(ProcessList)
% 1: add front
% 2: add back
% 3: detach front
% 4: detach back

out = zeros(length(ProcessList),1);

k_addpi = 1;
k_detpi = 1;
k_addalkyl = 1;
k_detalkyl = 1;
k_addpialkyl = 1;
k_detpialkyl = 1;

out(1) = k_addpi*ProcessList(1,2);
out(2) = k_detpi*ProcessList(2,2);
out(3) = k_addalkyl*ProcessList(3,2);
out(4) = k_detalkyl*ProcessList(4,2);
out(5) = k_addpialkyl*ProcessList(5,2);
out(6) = k_detpialkyl*ProcessList(6,2);

end

% function out = add_pi()
% k = 1;
% out = k;
% end
% 
% function out = det_pi()
% k = 1;
% out = k;
% end
% 
% function out = det_front_rate(Stack)
% global T kb
% d0 = intrinsic_detach();
% [overlap,OV] = find_overlap(Stack,length(Stack),length(Stack)-1);
% Eov = overlap_E();
% out = d0*exp(-overlap*Eov/(kb*T));
% end
% 
% function out = det_back_rate(Stack)
% global T kb
% d0 = intrinsic_detach();
% [overlap,OV] = find_overlap(Stack,1,2);
% Eov = overlap_E();
% out = d0*exp(-overlap*Eov/(kb*T));
% end
% 
% function out = intrinsic_detach()
% out = 8E11;
% end
% 
% function out = overlap_E()
% out = 1.23E-21;
% end


function [overlap_length,overlap_vector] = find_overlap(i,j,layer)
global Stack
overlap_vector = [max(Stack(i,1,layer),Stack(j,1,layer)), min(Stack(i,2,layer),Stack(j,2,layer))]; % find overlap of previous two chains
overlap_length = overlap_vector(2)-overlap_vector(1)+1;

end

% function DPs = redefine_add(new_chain)
% global DPdist1
% x = find(DPdist1==new_chain(2)-new_chain(1)+1);
% DPdist1(x(discretesample(x,1))) = [];
% DPs = tabulate(DPdist1);
% 
% end

%% Distribution Changes

function DPs = redefine_det(det_chain)

global DPdist1
DPdist1 = [DPdist1; det_chain(2)-det_chain(1)+1];

end

function ChainLength = DPSample()

global DPdist1
NumFreeChains = size(DPdist1,1);    % How many chains are solvated, total
Which = randi(NumFreeChains,1);     % Pick one
ChainLength = DPdist1(Which,1);     % Find its length
DPdist1(Which,:) = [];              % Remove it from the distribution

end
