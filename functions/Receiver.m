classdef Receiver < handle
    %{
        -------------------------------------------------------------------
        Microphone operations:
        -------------------------------------------------------------------
        __"addMic" [x y z] [theta phi r]
        adds a microphone, with certain coordinates and a direction vector
        
        __"rmMic" #mic
        removes microphone in position #mic
    
        __"rmAllMic"
        removes microphone in position #mic

        __"cMicPos" #mic [x y z]
        changes the position coordinates of the #mic microphone

        __"cMicDir" #mic [theta phi r]
        changes the direction vector of the #mic microphone

        __"cMic" #mic [x y z] [theta phi r]
        changes the direction vector of the #mic microphone
        changes the direction vector of the #mic microphone

        __"tf" @(TF callback)
        changes all microphones transfer function
        
        -------------------------------------------------------------------
        Walls/Surfaces operations:
        -------------------------------------------------------------------
        __"addWall" [x1 y1 z1] [x2 y2 z2] [x3 y3 z3] [x4 y4 z4] RC
        adds a wall/surface with the coordinates of 4 points, indicating
        its reflection coefficient RC

        __"rmWall" #wall
        removes wall in position #wall <=> removes a total of 4 consecutive
        points and 1 reflection coefficient
    
        __"rmAllWalls"
        removes wall in position #wall <=> removes a total of 4 consecutive
        points and 1 reflection coefficient

        __"cWallPos" #wall [x1 y1 z1] [x2 y2 z2] [x3 y3 z3] [x4 y4 z4]
        changes the 4 point coordinates of #wall wall/surface

        __"cWallRC" #wall RC
        changes the reflection coefficient of #wall wall/surface
    
        __"cWall" #wall [x1 y1 z1] [x2 y2 z2] [x3 y3 z3] [x4 y4 z4] RC
        changes the 4 point coordinates of #wall wall/surface
        changes the reflection coefficient of #wall wall/surface
    %}
    
    % Public, tunable properties
    properties
        mic
        s
    end

    % Pre-computed constants
    properties(Access = private)
    end

    methods
        % default constructor
        function obj = Receiver(struc, varargin)
            % using 'recstruct.mat':
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % struc = *.mat file that contains the structure.
            % The structure must have the following field tree:
            % ReceiverName  ___ mic__. pos
            %                     |__. dir
            %                     |__. tf
            %               \__ s __. pos
            %                    |__. rc
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % [1] struc = receiver structure (saved file)
            load(struc);
            obj.mic = rec.mic;
            obj.s = rec.s;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % [2] dmf = distance multiplying factor
            if size(varargin,2) == 1
                dmf = varargin{1};
                % microphone
                N = size(obj.mic.pos,1);
                for i=1:N
                    obj.mic.pos(i,2) = obj.mic.pos(i,2)*dmf;
                end
                % head
                
            elseif size(varargin,2) > 1
                error('MyComponent:incorrectType',...
                    'Receiver(structure)\nor Receiver(structure,dmf)');
            end
        end
        
        %____________________________Microphone____________________________
        function obj = addMic(obj, coord, dir)
            % structure check
            if (length(coord) ~= 3)
                error('Coordinates must respect the format [x y z].');
            elseif (~isnumeric(coord))
                error('Coordinates must contain numeric values ONLY.');
            elseif (length(dir) ~= 2)
                error('Direction vector must respect the format [theta phi].');
            elseif (~isnumeric(dir))
                error('Direction vector must contain numeric values ONLY.');
            end
            
            % values integrity
            if (abs(dir(1)) > pi)
                error('Theta must be between [-pi pi] radians.');
            elseif (abs(dir(2)) > pi/2)
                error('Phi must be between [-pi/2 pi/2] radians.');
            end
            
            % ops
            N = size(obj.mic.pos,1);
            obj.mic.pos(N+1,:) = coord;
        end
        
        function obj = rmMic(obj, index)
            % structure check
            N = size(obj.mic.pos,1);
            if (N == 0)
                error('No microphone to remove.');
            elseif ((length(index)~=1) && (~isnumeric(index)))
                error('Index must be a numeric scalar value.');
            end
            
            % integrity check
            if(index <= 0)
                error('Index must be a positive value (this is MATLAB).');
            elseif(index > N)
                error('Index should be smaller or equal to %d.', N);
            end
            
            % ops
            obj.mic.pos(index,:) = [];
            obj.mic.dir(index,:) = [];
        end
        
        function obj = rmAllMic(obj)
            obj.mic.pos = [];
            obj.mic.dir = [];
        end
        
        function obj = cMicPos(obj, index, pos)
        %__"cMicPos" #mic [x y z]
        %changes the position coordinates of the #mic microphone
            checkIndex(obj,index);
            checkCoordinates(obj,pos);
            obj.mic.pos(index,:) = pos;
        end
        
        function obj = cMicDir(obj, index, dir)
            % structure check
            if ((length(index)~=1) && (~isnumeric(index)))
                error('Index must be a numeric scalar value.');
            elseif (length(dir) ~= 2)
                error('Direction vector must respect the format [theta phi].');
            elseif (~isnumeric(dir))
                error('Director vector must contain numeric values ONLY.');
            end
            
            % values integrity
            % integrity check
            if(index <= 0)
                error('Index must be a positive value (this is MATLAB).');
            elseif(index > N)
                error('Index should be smaller or equal to %d.', N);
            elseif(abs(dir(1)) > pi)
                error('Theta must be between [-pi pi] radians.');
            elseif (abs(dir(2)) > pi/2)
                error('Phi must be between [-pi/2 pi/2] radians.');
            end
            
            obj.mic.dir(index,:) = dir;
        end
        
        function obj = cMic(obj, index, pos, dir)
            obj.cMicPos(index, pos);
            obj.cMicDir(index, dir);
        end
        
        function obj = cMicTF(obj,callback)
            obj.mic.tf = callback;
        end
        
        %__________________________Walls/Surfaces__________________________
        function obj = addWalls(obj, varargin)
            % structure check
            if (length(varargin) ~= 5)
                error('Input = [x1 y1 z1] [x2 y2 z2] [x3 y3 z3] [x4 y4 z4] RC');
            end
            
            for i=1:4
                if (length(varargin{i}) ~= 3)
                    error('Coordinates must respect the format [x y z].');
                elseif (~isnumeric(varargin{i}))
                    error('Coordinates must contain numeric values ONLY.');
                end
            end
            
            if (length(varargin{5}) ~= 1)
                error('Reflection Coefficient must be a scalar value.');
            end
            
            % ops
            N = size(obj.s.pos,1);
            for i=1:4
                obj.s.pos(N+i,:) = varargin{i};
            end
            obj.s.rc = [obj.s.rc; varargin{5}];
        end
        
        function obj = rmWall(obj, index)
           % structure check
            N = size(obj.s.pos,1);
            if (N < 4)
                error('No walls to remove.');
            elseif ((length(index)~=1) && (~isnumeric(index)))
                error('Index must be a numeric scalar value.');
            end
            
            % integrity check
            if(index <= 0)
                error('Index must be a positive value (this is MATLAB).');
            elseif(index > floor(N/4))
                error('Index should be smaller or equal to %d.', floor(N/4));
            end
            
            % ops
            obj.s.pos(index*4-3:index*4,:) = [];
            obj.s.rc(index,:) = []; 
        end
        
        function obj = rmAllWalls(obj)
            obj.s.pos = [];
            obj.s.rc = []; 
        end
        
        function  obj = cWallPos(obj, varargin)
           % structure check
            if (length(varargin) ~= 5)
                error('Input = #wall [x1 y1 z1] [x2 y2 z2] [x3 y3 z3] [x4 y4 z4]');
            end
            
            N = size(obj.s.pos,1);
            index = varargin{1};
            if (N < 4)
                error('No walls to change.');
            elseif ((length(index)~=1) && (~isnumeric(index)))
                error('Index must be a numeric scalar value.');
            end
            
            for i=1:4
                if (length(varargin{i+1}) ~= 3)
                    error('Coordinates must respect the format [x y z].');
                elseif (~isnumeric(varargin{i+1}))
                    error('Coordinates must contain numeric values ONLY.');
                end
            end
            
            % integrity check
            if(index <= 0)
                error('Index must be a positive value (this is MATLAB).');
            elseif(index > floor(N/4))
                error('Index should be smaller or equal to %d.', floor(N/4));
            end
            
            % ops
            for i=1:4
                obj.s.pos(index*4-4+i,:) = varargin{i+1};
            end
        end
        
        function obj = cWallRC(obj, index, RC)
            N = size(obj.s.pos,1);
            if (N < 4)
                error('No walls to change the RC.');
            elseif ((length(index)~=1) && (~isnumeric(index)))
                error('Index must be a numeric scalar value.');
            end
            
            % integrity check
            if(index <= 0)
                error('Index must be a positive value (this is MATLAB).');
            elseif(index > floor(N/4))
                error('Index should be smaller or equal to %d.', floor(N/4));
            end
            
            obj.s.rc(index) = RC;
        end
        
        function obj = cWall(obj, varargin)
            index = varargin{1};
            s1 = varargin{2};
            s2 = varargin{3};
            s3 = varargin{4};
            s4 = varargin{5};
            obj.cWallPos(index, s1, s2, s3, s4);
            obj.cWallRC(index,varargin{6});
        end
        
        function checkIndex(obj,index)
            % structure check
            Ni = length(index); % index length
            if      Ni~=1
                    error('Index must be a scalar value.');
            elseif  ~isnumeric(index)
                    error('Index must be a numeric value.');
            end
            % integrity check
            Na = size(obj.mic.pos,1); % number of microphones arrays
            if      index <= 0
                error('Index must be a positive value (this is MATLAB).');
            elseif  index > Na
                error('Index should be smaller or equal to %d.', N);
            end
        end
        
        function checkCoordinates(obj,pos)
            % number of points in each microphone array
            Nc = size(obj.mic.pos,2); 
            % check
            if      Nc ~= 3
                error('Coordinates must respect the format [x y z].');
            elseif  ~isnumeric(pos)
                error('Coordinates must contain numeric values ONLY.');
            end
        end
        
    end
    
    methods (Static)
        
    end
end
