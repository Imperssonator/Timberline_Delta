function [out, L] = Kinetic_Stack(DPdist)

clc
%% Stackemup
% This function takes a degree of polymerization distribution as an input, in the form of an Nx1
% vector:
% The output will be a Nx2xM matrix where N is the number of chains in the
% first (longest) stack and M is the number of alkyl layers in the stck

%% Hard Coded Parameters
global T kb DPdist1 Stack
DPdist1 = DPdist;
T = 298;
kb = 1.38E-23;
Iterations = 5500;
figure(4)
hist(DPdist1,20)    % Histogram of initial distribution
Stack = [0 DPSample()-1]; % stack is zero-indexed at the first monomer of the first chain to be picked

L = zeros(Iterations,20);
iter = 0;

while iter<Iterations 
    iter = iter+1;
    if isempty(DPdist1)
        break
    end
    if size(Stack,1)>1
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
        for i = 1:size(Stack,3)
            w = Stack(:,:,i);   
            L(iter,i) = length(w(~all(w==0,2),:));  % Number of non-zero rows in the matrix w ie length of that layer of the stack
        end
        disp(size(Stack,1))
    else
        Stack = initiate();
    end
    
       
%     if mod(iter,10000) == 0
%         figure(1)
%         hold on
%         stackplot(Stack,iter)
%         drawnow
        
       
end
out = Stack;
L(:,size(Stack,3):end) = [];

figure(1)
stackplot(Stack)
hold off

StartDist = tabulate(DPdist);
First = StartDist(1,1); Last = StartDist(end,1);
EndDist = tabulate(DPdist1);
InnerFirst = EndDist(1,1); InnerLast = EndDist(end,1);
Prefix = (First:InnerFirst-1)';
Suffix = (InnerLast+1:Last)';
Prefix = [Prefix, zeros(size(Prefix,1),2)];
Suffix = [Suffix, zeros(size(Suffix,1),2)];
EndDist = [Prefix; EndDist; Suffix];                % Add zeros for the lengths that were totally depleted in the simulation

UsedChains = [StartDist(:,1), StartDist(:,2:3)-EndDist(:,2:3)];
% figure
% bar(UsedChains(:,1),UsedChains(:,2))
% 
figure(2)
plot(L)

figure (3)
hist(DPdist1,20);           % Histogram of remaining free chains
hold off
 
end

function out = initiate()
%% Initiate
% Add a chain if there is only 1
global Stack
DP = DPSample();
OV = Stack(1,:,1);

hitstack = randi(OV(2)-OV(1)+1)+OV(1)-1; % where does new chain hit stack
hitchain = randi(DP);                    % which monomer of new chain hits stack
new_chain = [hitstack-hitchain+1, hitstack+DP-hitchain];
Stack(2,:,1) = new_chain;
out = Stack;
end

function ChainLength = DPSample()
% Picks a random chain from the distribution
global DPdist1
NumFreeChains = size(DPdist1,1);    % How many chains are solvated, total
Which = randi(NumFreeChains,1);     % Pick one
ChainLength = DPdist1(Which,1);     % Find its length
DPdist1(Which,:) = [];              % Remove it from the distribution

end

function out = find_process()
% This function gives a list of the number of ways in which each process
% can take place. The output is a 6x2 matrix
global Stack occupied_top occupied_bottom growth_sites removables FirstLayer
[m, ~, p] = size(Stack);
out = [1 0; 2 0; 3 0; 4 0; 5 0; 6 0];
out(1,2) = 2;   % Two pi chains can be added in either direction

FirstLayer = locate_stack();    % Locates the first layer in the stack (not the top layer, the FIRST layer that ever started to grow)
if size(Stack,3) == 1
    out(2,2) = 2;
else
    if FirstLayer ~= 1
        if FirstLayer ~= size(Stack,3)
            if Stack(1,:,FirstLayer-1) == [0 0] & Stack(1,:,FirstLayer+1) == [0 0]
                out(2,2) = out(2,2)+1;
            end
            if Stack(size(Stack,1),:,FirstLayer-1) == [0 0] & Stack(size(Stack,1),:,FirstLayer+1) == [0 0]
                out(2,2) = out(2,2)+1;
            end
        else
            if Stack(1,:,FirstLayer-1) == [0 0]
                out(2,2) = out(2,2)+1;
            end
            if Stack(size(Stack,1),:,FirstLayer-1) == [0 0] 
                out(2,2) = out(2,2)+1;
            end
        end
    else
        if Stack(1,:,FirstLayer+1) == [0 0]
            out(2,2) = out(2,2)+1;
        end
        if Stack(size(Stack,1),:,FirstLayer+1) == [0 0]
            out(2,2) = out(2,2)+1;
        end
    end
end
top_layer = Stack(:,:,1);   % Top layer of the stack
occupied_top = find_sites(top_layer);   % Number of sites occupied in the top 
% layer to get the number of possible sites a new alkyl layer can nucleate
if isempty(occupied_top) == 1
    occupied_top = find_sites(Stack(:,:,2));
    Stack(:,:,1) = [];
end
    
bottom_layer = Stack(:,:,end);  % Bottom layer of the stack
occupied_bottom = find_sites(bottom_layer);
if isempty(occupied_bottom) == 1
    occupied_bottom = find_sites(Stack(:,:,end-1));
    Stack(:,:,end) = [];
end
out(3,2) = length(occupied_top) + length(occupied_bottom);  % Number of places where new alkyl chain can begin

if length(occupied_top) == 1 & length(occupied_bottom) == 1 & p>1
    out(4,2) = 2;   % Process 4 is possible only when there is just one chain in the alkyl layer
elseif length(occupied_top) == 1 || length(occupied_bottom) == 1 & p>1
    out(4,2) = 1;
end
% Number of ways in which single alkyl stack chain can be removed

growth_sites = [];
removables = [];
for i = 1:FirstLayer-1
    occupied = find_sites(Stack(:,:,i));    % Indices of the occupied sites in the ith layer
    if isempty(occupied) == 1
        continue
    end
    if occupied(1) ~= 1 & Stack(occupied(1)-1,:,i+1) ~= [0 0]  % the first occupied site isn't the corner most and there is a chain in the layer below to allow for growth
        out(5,2) = out(5,2) + 1;
        growth_sites = [growth_sites; [occupied(1) i 1]];   % Keep track of possible growth site
    end
    if occupied(end) ~= m
        if Stack(occupied(end)+1,:,i+1) ~= [0 0]  %Same as above, only checking the other corner
            out(5,2) = out(5,2) + 1;
            growth_sites = [growth_sites; [occupied(end) i 2]];
        end
    end
    if i == 1 & size(occupied)>1   % if there is more than one chain in the top-most the layer, two corners can possibly be removed
        out(6,2) = out(6,2)+2;
        removables = [removables; [occupied(1),i]; [occupied(end),i]];  % Keep track of possible "removables"
    end
    if i>1 & Stack(occupied(end),:,i-1) == [0 0]   % Corners of the other layers not bound by more than two nearest chains
        out(6,2) = out(6,2) + 1;
        removables = [removables; [occupied(end),i]];
    end
    if i>1 & Stack(occupied(1),:,i-1) == [0 0]
        out(6,2) = out(6,2) + 1;
        removables = [removables; [occupied(1),i]];
    end
end
% Now the same stuff for the layers below the first layer
for i = size(Stack,3):-1:FirstLayer+1
    occupied = find_sites(Stack(:,:,i));
    if isempty(occupied) == 1
        continue
    end
    if occupied(1) ~= 1 & Stack(occupied(1)-1,:,size(Stack,3)-1) ~= [0 0]
        out(5,2) = out(5,2) + 1;
        growth_sites = [growth_sites; [occupied(1) i 1]];
    end
    if occupied(end) ~= m & Stack(occupied(end)+1,:,size(Stack,3)-1) ~= [0 0]
        out(5,2) = out(5,2) + 1;
        growth_sites = [growth_sites; [occupied(end) i 2]];
    end
    if i == size(Stack,3) & size(occupied) > 1
        out(6,2) = out(6,2)+2;
        removables = [removables; [occupied(1) i]; [occupied(end),i]];
    end
    if i<size(Stack,3) & Stack(occupied(end),:,i+1) == [0 0]
        out(6,2) = out(6,2) + 1;
        removables = [removables; [occupied(end),i]];
    end
    if i<size(Stack,3) & Stack(occupied(1),:,i+1) == [0 0]
        out(6,2) = out(6,2) + 1;
        removables = [removables; [occupied(1),i]];
    end
        
end
% Number of ways in which pi-alkyl chain can be added or removed
end

function out = find_sites(single_layer)
out = [];
find_zero = single_layer==0;    % Assigns zeros in the place of non zero elements and 1s in the place of zeros
find_nonzero = ~all(find_zero,2);   % Vector with 1s where the corresponding row in find_zero has at least one zero and 0 otherwise
out = find(find_nonzero);   % Indices of rows which have at least one non zero element ie. occupied sites
end

function out = locate_stack()   % This function finds the first layer that ever started to grow by finding the layer with zero unoccupied sites
global Stack
out = [];
for j = 1:size(Stack,3)
    occupied = find_sites(Stack(:,:,j));
    if length(occupied) == size(Stack,1)
        out = j;
        break
    end
end
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
% Output is the new stack
global occupied_top occupied_bottom growth_sites 
q = randi(List(p,2));
if p == 1   % add free chain to pi-stack, can be done only to first layer
    if q == 1   % choose between front and back
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
    if q<=length(occupied_top)
        out = add_top(q);
    else
        out = add_bottom(q-length(occupied_top));
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

function out = add_front()  % Adds a new chain to the front of the stack. 
                            % Output is the new stack
global Stack FirstLayer
m = size(Stack,1);
new_chain = collide(FirstLayer,size(Stack,1),size(Stack,1)-1);
Stack(m+1,:,:) = 0;
Stack(m+1,:,FirstLayer) = new_chain;
out = Stack;

end

function out = add_back()
global Stack FirstLayer
p = size(Stack,3);
new_chain = collide(FirstLayer,1,2);
Stack = cat(1,zeros(1,2,p),Stack);
Stack(1,:,FirstLayer) = new_chain;
out = Stack;

end

function out = collide(layer,front,support)
%% collide
% takes a chain from DPs and collides a random
% monomer with the frontier chain (front, 1x2 vector) and any location
% where front overlaps with support (1x2 vector)

% the output is a 1x2 vector
% pi stacking only

DP = DPSample();
[~,OV] = find_overlap(front,support,layer,layer); % find overlap of previous two chains
hitstack = randi((OV(2)-OV(1)+1))+OV(1)-1; % where does new chain hit stack
hitchain = randi(DP);                    % which monomer of new chain hits stack

out = [hitstack-hitchain+1, hitstack+DP-hitchain];

end

function [overlap_length,overlap_vector] = find_overlap(i,j,layer1,layer2)
global Stack
overlap_vector = [max(Stack(i,1,layer1),Stack(j,1,layer2)), min(Stack(i,2,layer1),Stack(j,2,layer2))]; % find overlap of previous two chains
overlap_length = overlap_vector(2)-overlap_vector(1)+1;
end

function out = det_front()
global Stack FirstLayer
redefine_det(Stack(end,:,FirstLayer));
Stack(end,:,:) = [];
out = Stack;

end

function out = det_back()
global Stack FirstLayer
redefine_det(Stack(1,:,FirstLayer));
Stack(1,:,:) = [];
out = Stack;

end


% function out = pick_pol()
% global DPs
% out = DPs(discretesample(DPs(:,3),1),1);
% end



function out = add_top(q)
global Stack occupied_top
Stack = cat(3, collide_alkyl(occupied_top(q),1),Stack);
out = Stack;
end

function out = add_bottom(q)
global Stack occupied_bottom
Stack = cat(3, Stack, collide_alkyl(occupied_bottom(q),size(Stack,3)));
out  = Stack;
end

function out = collide_alkyl(chain,layer)   % Begins a new alkyl layer with a chain attached to position mentioned
global Stack
[m, n, ~] = size(Stack);
out = zeros(m,n);
DP = DPSample();
hitstack = randi(Stack(chain,2,layer) - Stack(chain,1,layer) + 1) + Stack(chain,1,layer) - 1;
hitchain = randi(DP);

new_chain = [hitstack-hitchain+1, hitstack+DP-hitchain];
out(chain,:) = new_chain;
end

function out = det_top()
global Stack occupied_top
redefine_det(Stack(occupied_top(1),:,1));
Stack(:,:,1) = [];
out = Stack;
end

function out = det_bottom()
global Stack occupied_bottom
redefine_det(Stack(occupied_bottom(1),:,end));
Stack(:,:,end) = [];
out = Stack;
end

function out = add_pialkyl(front_chain)
global Stack FirstLayer
new_chain = collide_pialkyl(front_chain);
if front_chain(2)< FirstLayer
     SupportLayer = front_chain(2)+1;
else
    SupportLayer = front_chain(2)-1;
end
if front_chain(3) == 1
    Stack(front_chain(1)-1,:,front_chain(2)) = new_chain;
    [OL,OV] = find_overlap(front_chain(1)-1,front_chain(1)-1,front_chain(2),SupportLayer);
else
    Stack(front_chain(1)+1,:,front_chain(2)) = new_chain;
    [OL,OV] = find_overlap(front_chain(1)+1,front_chain(1)+1,front_chain(2),SupportLayer);
end
if OL<=0
    Stack = add_pialkyl(front_chain);
end
out = Stack;
end

function out = collide_pialkyl(front_chain)
global Stack 
DP = DPSample();
hitstack = randi(Stack(front_chain(1),2,front_chain(2))-Stack(front_chain(1),1,front_chain(2))+1) + Stack(front_chain(1),1,front_chain(2)) - 1;
hitchain = randi(DP);

out = [hitstack-hitchain+1, hitstack+DP-hitchain];

end

function out = det_pialkyl(q)
global Stack removables 
det = removables(q,:);
redefine_det(Stack(det(1),:,det(2)));
Stack(det(1),:,det(2)) = [0 0];
out =  Stack;
end

%% RATES
% Everything involved with getting rates is in this section
function out = get_rates(ProcessList)

out = zeros(length(ProcessList),1);

% Rate constants for each of the processes
k_addpi = 1.2;
k_detpi = 1;
k_addalkyl = 0.0001;
k_detalkyl = 2;
k_addpialkyl = 1.5;
k_detpialkyl = 0.8;

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

