function CommanderMissile
 
    globals.explodeSpeed = 0.03;
    globals.explodeRadius = 0.1;
    globals.yPixels = 360;
    globals.xPixels = 640;
    globals.locations = [[.87 .9 .92 .92 .86 .92 .9 .92 .87]' (((0:8)+0.5)/9)']';
    globals.ax  = gca;
    globals.interceptorSpeed = 0.5;
    globals.smartTurnRadius = 0.05;
    
    im = [];
    controlAxis = [];
    
    state = 'I';
    interrupt = 0;
    tic;
    endTime = inf;
    for iLoop = 1:1e12
                
        time = toc;
                   
        if (~ishandle(globals.ax))
            break;
        end
        
        if (time>endTime)
            state = 'S';
        end
        
        switch state
            
            case 'I'
                
                im = image(uint8(255 * ones(globals.yPixels,globals.xPixels,3)));
                imAxis = gca;
                set(gcf,'pointer','crosshair');
                if (isempty(controlAxis))
                    controlAxis = axes;
                    controlAxis.Color = 'none';
                    controlAxis.XLim = globals.ax.XLim;
                    controlAxis.YLim = globals.ax.YLim;
                    controlAxis.YDir = 'reverse';
                    controlAxis.ButtonDownFcn = @AxisCallback;
                    globals.controlAxis = controlAxis;
                end
                
                imAxis.Position = [0 0 1 1];
                imAxis.XColor = 'none';
                imAxis.YColor = 'none';
                controlAxis.Position = [0 0 1 1];
                controlAxis.XColor = 'none';
                controlAxis.YColor = 'none';
                
                globals.yMissile = LaunchNoise;
                globals.fMissile = 8192;
                
                globals.yExplosion = ExplosionNoise;
                globals.fExplosion = 8192;
                
                globals.ySiren = SirenNoise;
                globals.fSiren = 8192;
                
                HandleCities([],[],inf,[],globals);
                [threatSets,colorSchemes] = MakeWaves;
                wave = 1;
                
                state = 'M';
                
            case 'M'

                % 'missile bomber satellite'
                threatSet = threatSets{wave};
                [missileArray, bomberArray] = MakeThreats(threatSet,globals);

                colorScheme = colorSchemes{wave};
                background = MakeBackground(globals,colorScheme);
                playImage = background;
                
                state = 'P';
                
            case 'P'
                
                explodeArray = [];
                [playImage,explodeArray] = HandleMissiles(playImage,explodeArray,-1,[],background,globals);
                [playImage,explodeArray] = HandleBombers(playImage,explodeArray,-1,[],background,globals);
                [playImage,explodeArray] = HandleSmartBombs(playImage,explodeArray,-1,[],background,globals);
                [playImage,explodeArray] = HandleExplosions(playImage,explodeArray,-1,background,globals);
                [playImage,explodeArray] = HandleInterceptor(playImage,explodeArray,-1,background,globals);
                [playImage,explodeArray] = HandleCities(playImage,explodeArray,-1,background,globals);
                
                for idx = 1:size(missileArray,2)
                    [playImage,explodeArray] = HandleMissiles(playImage,explodeArray,inf,missileArray(:,idx),background,globals);
                end
                for idx = 1:size(bomberArray,2)
                    [playImage,explodeArray] = HandleBombers(playImage,explodeArray,inf,bomberArray(:,idx),background,globals);
                end
                soundsc(globals.ySiren,globals.fSiren);
                state = 'W';

            case 'W'
                
                figure(gcf);
                im.CData = playImage;
                drawnow;
                   
                % Wait before start of wave
                if (time < 3)
                    pause(0.1); 
                else
                    tic
                    state = 'R';
                end
                
            case 'R'
                
                switch mod(interrupt,8)
                    case 0
                        im.CData = playImage;
                        drawnow;
                    case 1
                        [playImage,explodeArray] = HandleMissiles(playImage,explodeArray,time,[],background,globals);
                    case 2
                        [playImage,explodeArray] = HandleExplosions(playImage,explodeArray,time,background,globals);
                    case 3
                        [playImage,explodeArray] = HandleInterceptor(playImage,explodeArray,time,background,globals);
                    case 4
                        [playImage,explodeArray] = HandleSmartBombs(playImage,explodeArray,time,[],background,globals);
                    case 5
                        [playImage,explodeArray] = HandleBombers(playImage,explodeArray,time,[],background,globals);
                    case 6  
                        [playImage,explodeArray] = HandleCities(playImage,explodeArray,time,background,globals);
                    case 7
                        numLeftMissilesLeft = HandleMissiles(playImage,explodeArray,nan,[],background,globals);
                        numLeftSmartBombsLeft = HandleSmartBombs(playImage,explodeArray,nan,[],background,globals);
                        numLeftBombersLeft = HandleBombers(playImage,explodeArray,nan,[],background,globals);
                        if ((numLeftMissilesLeft+numLeftSmartBombsLeft+numLeftBombersLeft)==0)
                            explodeTime = 2*globals.explodeRadius/globals.explodeSpeed;
                            if (endTime > time + explodeTime)
                                endTime = time + explodeTime;
                            end
                        else
                            endTime = inf;
                        end
                end
                interrupt = interrupt + 1;

            case 'S'
                
                tic;
                wave = wave + 1;
                if ( HandleCities(playImage,explodeArray,nan,background,globals)==0 || ...
                        wave>length(threatSets) )
                    break;
                end
                tic;
                state = 'M';
                
        end
        
    end
    
    disp('THE END');
    close all;

end

function background = MakeBackground(globals,colorScheme)

    background = zeros(globals.yPixels,globals.xPixels,3,'uint8');
    background(:,:,1) = colorScheme(1,1);
    background(:,:,2) = colorScheme(1,2);
    background(:,:,3) = colorScheme(1,3);
    
    dHeightA = 0.03;
    dHeightB = 0.02;
    dWidthA = 0.02;
    dWidthB = 0.03;
    ground = [
        globals.locations(1,2)+[0 0] ...
        globals.locations(1,1)+[-dHeightA -dHeightA] ...
        globals.locations(1,2)+[0 0 -dHeightB -dHeightB] ...
        globals.locations(1,3)+[0 0 -dHeightB -dHeightB] ...
        globals.locations(1,4)+[0 0 -dHeightB -dHeightB 0 0] ...
        globals.locations(1,5)+[-dHeightA -dHeightA] ...
        globals.locations(1,6)+[0 0 -dHeightB -dHeightB] ...
        globals.locations(1,7)+[0 0 -dHeightB -dHeightB] ...
        globals.locations(1,8)+[0 0 -dHeightB -dHeightB 0 0] ...
        globals.locations(1,9)+[-dHeightA -dHeightA] ...
        globals.locations(1,8)+[0 0];
        -1 globals.locations(2,1)+[-dWidthA-dWidthB -dWidthB dWidthB dWidthB+dWidthA] ...
        globals.locations(2,2)+[-dWidthA-dWidthB -dWidthB dWidthB dWidthB+dWidthA] ...
        globals.locations(2,3)+[-dWidthA-dWidthB -dWidthB dWidthB dWidthB+dWidthA] ...
        globals.locations(2,4)+[-dWidthA-dWidthB -dWidthB dWidthB dWidthB+dWidthA] ...
        globals.locations(2,5)+[-dWidthA-dWidthB -dWidthB dWidthB dWidthB+dWidthA] ...
        globals.locations(2,6)+[-dWidthA-dWidthB -dWidthB dWidthB dWidthB+dWidthA] ...
        globals.locations(2,7)+[-dWidthA-dWidthB -dWidthB dWidthB dWidthB+dWidthA] ...
        globals.locations(2,8)+[-dWidthA-dWidthB -dWidthB dWidthB dWidthB+dWidthA] ...
        globals.locations(2,9)+[-dWidthA-dWidthB -dWidthB dWidthB dWidthB+dWidthA] 2];
    ground(1,:) = 1 + ground(1,:)*(globals.yPixels-1);
    ground(2,:) = 1 + ground(2,:)*(globals.xPixels-1);
    groundInterpolator = griddedInterpolant(ground(2,:),ground(1,:));
    yPixels = 1:globals.yPixels;
    for xPixel = 1:globals.xPixels
        yPixelHeight = groundInterpolator(xPixel);
        idx = yPixels>=yPixelHeight;
        background = SubColor(background,yPixels(idx)',xPixel*ones(sum(idx),1),colorScheme(2,:));
    end

end

function [missileArray,bomberArray] = MakeThreats(threatSet,globals)

    missileArray = [];
    bomberArray = [];
    
    targetLocations = [
        1 + globals.locations(1,:)*(globals.yPixels-1);
        1 + globals.locations(2,:)*(globals.xPixels-1)];
            
    for iThreat = 1:length(threatSet)
        
        threatCommands = split(threatSet{iThreat});
        if (any(contains(threatCommands,'missile')))
            idxMissiles = find(contains(threatCommands,'missile'))+1;
            idxTime = find(contains(threatCommands,'time'))+1;
            idxSpeed = find(contains(threatCommands,'speed'))+1;
            numWarheads = str2double(threatCommands{idxMissiles(1)});
            startTime = str2double(threatCommands{idxTime(1)});
            speed = globals.yPixels*str2double(threatCommands{idxSpeed(1)});
            startPosition = [0 (globals.xPixels-1)*rand(1)]';
            endPosition = targetLocations(:,ceil(9*rand(1))) + globals.yPixels*0.01*randn(2,1);
            distance = sqrt(sum((endPosition-startPosition).^2));
            endTime = startTime + distance/speed;
            missileArray(:,end+1) = [startTime;endTime;startPosition;endPosition;numWarheads];
        end
        if (any(contains(threatCommands,'bomber')))
            idxBombers = find(contains(threatCommands,'bomber'))+1;
            idxTime = find(contains(threatCommands,'time'))+1;
            idxSpeed = find(contains(threatCommands,'speed'))+1;
            numBombs = str2double(threatCommands{idxBombers(1)});
            startTime = str2double(threatCommands{idxTime(1)});
            speed = globals.yPixels*str2double(threatCommands{idxSpeed(1)});
            flyingLeft = rand(1)>0.5;
            if (flyingLeft)
                startYx = [globals.yPixels*(0.25*rand(1)+0.25) globals.xPixels]';
                stopYx = [startYx(1) 1]';
            else
                startYx = [globals.yPixels*(0.25*rand(1)+0.25) 1]';
                stopYx = [startYx(1) globals.xPixels]';
            end
            distance = globals.xPixels-1;
            endTime = startTime + distance/speed;
            bomberArray(:,end+1) = [startTime;endTime;startYx;stopYx;numBombs];
        end
    end
    
end

function [playImage,explodeArray] = HandleMissiles(playImage,explodeArray,time,missleArray,background,globals)

    persistent missileStructs lastTime;

    if ( time<0 )
        
        missileStructs = [];
        
    elseif ( isinf(time) )
        
        yPixels = globals.yPixels;
        xPixels = globals.xPixels;
        startTime = missleArray(1);
        stopTime = missleArray(2);
        startYx = missleArray(3:4);
        stopYx = missleArray(5:6);
        [pixelYx,pixelTimes] = LineToPixels(yPixels,xPixels,startTime,stopTime,startYx,stopYx);
        missileStruct.pixelYx = pixelYx;
        missileStruct.pixelTimes = pixelTimes;
        missileStruct.yInterpolator = ...
            griddedInterpolant([startTime stopTime],[startYx(1) stopYx(1)]);
        missileStruct.xInterpolator = ...
            griddedInterpolant([startTime stopTime],[startYx(2) stopYx(2)]);
        missileStruct.numMirvs = missleArray(7);
        missileStruct.timeMirvs = startTime + (stopTime-startTime)*(0.25 + 0.25*rand(1));
        if ( isempty(missileStructs) )
            missileStructs = missileStruct;
        else
            missileStructs(end+1) = missileStruct;
        end
                        
        lastTime = time;
        
    elseif (isnan(time))
        
        playImage = length(missileStructs);
        
    else
        
        missileArray = [];
        goBoomArray = zeros(length(missileStructs),2);
        for idxMissile = 1:length(missileStructs)
            
            posYx = [missileStructs(idxMissile).yInterpolator(time) missileStructs(idxMissile).xInterpolator(time)]';
            goBoom = NaN(1,2);
            
            if ( time>=missileStructs(idxMissile).yInterpolator.GridVectors{1}(end) )
                goBoom = [missileStructs(idxMissile).yInterpolator.Values(end) missileStructs(idxMissile).xInterpolator.Values(end)];
            end
            
            if ( ~isempty(explodeArray) )
                distances = sqrt(sum((explodeArray(2:3,:)-posYx).^2,1));
                if (any((distances-explodeArray(4,:))<=0))
                    goBoom = posYx;
                end
            end
            
            goBoomArray(idxMissile,:) = goBoom;
 
            pixelYx = missileStructs(idxMissile).pixelYx;
            pixelTimes = missileStructs(idxMissile).pixelTimes;
            if ( isnan(goBoom(1)) )
                
                idx = pixelTimes>lastTime & pixelTimes<=time;                
                if ( any(idx) )
                    playImage = AddColor(playImage,pixelYx(idx,1),pixelYx(idx,2),[63 63 63]);
                end
                
                if ( time > missileStructs(idxMissile).timeMirvs && ...
                        missileStructs(idxMissile).numMirvs > 1 )
                    numMirvs = missileStructs(idxMissile).numMirvs;
                    targetLocations = [
                            1 + globals.locations(1,:)*(globals.yPixels-1);
                            1 + globals.locations(2,:)*(globals.xPixels-1)];
                    missileSpeed = ...
                        norm([missileStructs(idxMissile).yInterpolator.Values(end)-missileStructs(idxMissile).yInterpolator.Values(1) ...
                        missileStructs(idxMissile).xInterpolator.Values(end)-missileStructs(idxMissile).xInterpolator.Values(1)])/ ...
                        (missileStructs(idxMissile).yInterpolator.GridVectors{1}(end)-missileStructs(idxMissile).yInterpolator.GridVectors{1}(1));
                    for idx = 1:(numMirvs-1)
                        targetPosition = targetLocations(:,ceil(9*rand(1))) + globals.yPixels*0.01*randn(2,1);
                        distance = norm(targetPosition - posYx);
                        endTime = time + distance/missileSpeed;
                        missileArray(:,end+1) = [time;endTime;posYx;targetPosition;1];
                    end
                    missileStructs(idxMissile).numMirvs = 1;
                end
                
            else
                
                idx = pixelTimes<=time;                
                if ( any(idx) )
                    playImage = EraseColor(playImage,pixelYx(idx,1),pixelYx(idx,2),background);
                end
                
            end
            
        end

        idxBoom = find(~isnan(goBoomArray(:,1)))';
        for idxMissile = idxBoom
            if ( ~isnan(goBoomArray(idxMissile,1)) )
                explodeArray(:,end+1) = [time goBoomArray(idxMissile,1) goBoomArray(idxMissile,2) 0];
            end
        end
        missileStructs(idxBoom) = [];
        
        for idx = 1:size(missileArray,2)
            [playImage,explodeArray] = HandleMissiles(playImage,explodeArray,inf,missileArray(:,idx),background,globals);
        end
        
        lastTime = time;
        
    end

end

function [playImage,explodeArray] = HandleBombers(playImage,explodeArray,time,bomberArray,background,globals)

    persistent bomberStructs leftBomber rightBomber;

    if ( time<0 )
        
        bomberStructs = [];

        leftBomber = [0*(-11:9) -1+0*(-9:11) 1+0*(-10:7) -2+0*(8:12) -3+0*(10:14);
                        (-11:9)      (-9:11)     (-10:7)      (8:12)      (10:14)];
        leftBomber = horzcat(leftBomber, ...
            [2+0*(-4:1) 3+0*(-1:3) 4+0*(2:5) -2+0*(-1:3) -3+0*(2:5);
                 (-4:1)     (-1:3)     (2:5)      (-1:3)      (2:5)]);
        rightBomber = leftBomber;
        rightBomber(2,:) = -rightBomber(2,:);
            
    elseif ( isinf(time) )
        
            startTime = bomberArray(1);
            stopTime = bomberArray(2);
            startYx = bomberArray(3:4);
            stopYx = bomberArray(5:6);
            numBombs = bomberArray(7);
            bomberStruct.yInterpolator = griddedInterpolant([startTime stopTime],[startYx(1) stopYx(1)]);
            bomberStruct.xInterpolator = griddedInterpolant([startTime stopTime],[startYx(2) stopYx(2)]);
            bomberStruct.flyingLeft = startYx(2) > stopYx(2);
            flightTime = stopTime - startTime;
            if ( numBombs == 1 )
                bomberStruct.bombTimes = startTime + 0.5 * flightTime;
            elseif ( numBombs == 2 )
                bomberStruct.bombTimes = startTime + flightTime * [0.375 0.625];
            else
                bomberStruct.bombTimes = startTime + flightTime * (1/3:(1/3/(numBombs-1)):2/3);
            end
            bomberStruct.pixelsToDelete = [];
            bomberStruct.bomberSpeed = norm(stopYx - startYx)/flightTime;
            if ( bomberStruct.flyingLeft )
                bomberStruct.bomberSpeed = -bomberStruct.bomberSpeed;
            end
            if ( isempty(bomberStructs) )
                bomberStructs = bomberStruct;
            else
                bomberStructs(end+1) = bomberStruct;
            end
        
    elseif ( isnan(time) )
        
        playImage = length(bomberStructs);
        
    else

        goBoomArray = zeros(length(bomberStructs),2);
        bomberGoneArray = zeros(length(bomberStructs),1)>0;
        for idxBomber = 1:length(bomberStructs)

            if ( ~isempty(bomberStructs(idxBomber).pixelsToDelete) )
                playImage = EraseColor(playImage, ...
                    bomberStructs(idxBomber).pixelsToDelete(1,:),bomberStructs(idxBomber).pixelsToDelete(2,:),background);
            end
            
            posYx = [bomberStructs(idxBomber).yInterpolator(time);bomberStructs(idxBomber).xInterpolator(time)];
            goBoom = NaN(1,2);
            
            if ( ~isempty(explodeArray) )
                distances = sqrt(sum((explodeArray(2:3,:)-posYx).^2,1));
                if (any((distances-explodeArray(4,:))<=0))
                    goBoom = posYx;
                end
            end
            
            goBoomArray(idxBomber,:) = goBoom;
            
            if ( isnan(goBoom(1)) )
            
                if ( time>=bomberStructs(idxBomber).yInterpolator.GridVectors{1}(end) )
                    bomberGoneArray(idxBomber) = true;
                end
            
                if ( bomberStructs(idxBomber).flyingLeft )
                    yPixels = round(posYx(1)) + leftBomber(1,:);
                    xPixels = round(posYx(2)) + leftBomber(2,:);
                else
                    yPixels = round(posYx(1)) + rightBomber(1,:);
                    xPixels = round(posYx(2)) + rightBomber(2,:);
                end
                idx = yPixels>=1 & yPixels<=globals.yPixels & ...
                    xPixels>=1 & xPixels<=globals.xPixels;
                
                playImage = SubColor(playImage,yPixels(idx),xPixels(idx),[0 0 255]);
                bomberStructs(idxBomber).pixelsToDelete = [yPixels(idx);xPixels(idx)];
                
                if ( ~isempty(bomberStructs(idxBomber).bombTimes) )
                    bombTime = bomberStructs(idxBomber).bombTimes(1);
                    if ( time>=bombTime )
                        targetLocations = [
                            1 + globals.locations(1,:)*(globals.yPixels-1);
                            1 + globals.locations(2,:)*(globals.xPixels-1)];
                        targetPosition = targetLocations(:,ceil(9*rand(1))) + globals.yPixels*0.01*randn(2,1);
                        smartbombArray = [bombTime 0 bomberStructs(idxBomber).bomberSpeed ...
                            bomberStructs(idxBomber).yInterpolator(bombTime) bomberStructs(idxBomber).xInterpolator(bombTime) ...
                            targetPosition(1) targetPosition(2)]';
                        [playImage,explodeArray] = HandleSmartBombs(playImage,explodeArray,inf,smartbombArray,background,globals);
                        bomberStructs(idxBomber).bombTimes(1) = [];
                    end
                end
                
            end
 
        end
        
        idxBoom = find(~isnan(goBoomArray(:,1)))';
        for idxBomber = idxBoom
            if ( ~isnan(goBoomArray(idxBomber,1)) )
                explodeArray(:,end+1) = [time goBoomArray(idxBomber,1) goBoomArray(idxBomber,2) 0];
            end
        end
        bomberStructs(idxBoom) = [];
        
        if ( any(bomberGoneArray) )
            for idxBomber = find(bomberGoneArray)
                if ( ~isempty(bomberStructs(idxBomber).pixelsToDelete) )
                    playImage = EraseColor(playImage, ...
                        bomberStructs(idxBomber).pixelsToDelete(1,:),bomberStructs(idxBomber).pixelsToDelete(2,:),background);
                end
            end
        end
        bomberStructs(bomberGoneArray) = [];
        
    end

end

function [playImage,explodeArray] = HandleSmartBombs(playImage,explodeArray,time,smartbombArray,background,globals)

    persistent smartBombStructs pixelsToAdd;

    if ( time<0 )
        
        smartBombStructs = [];
        
        pixelsToAdd = [-2 -1 0 1 2  2  1 -1 -2;
            -2 -1 0 1 2 -2 -1  1  2];      
        
    elseif ( isinf(time) )
        
        velYx = smartbombArray(2:3);
        posYx = smartbombArray(4:5);
        targetYx =  smartbombArray(6:7);
        velVector = [velYx;0];
        posVector = [posYx;0];
        tgtVector = [targetYx;0];
        turnRadius = globals.yPixels*globals.smartTurnRadius;
        turnCenter = posVector + [turnRadius 0 0]';
        turnX = posVector - turnCenter;
        turnY = [turnRadius * velYx / norm(velYx);0];
        turnRate = norm(velVector)/turnRadius;
        
        smartBombStruct.posVector = posVector;
        smartBombStruct.tgtVector = tgtVector;
        smartBombStruct.pixelsToDelete = [];
        smartBombStruct.timeBomb = smartbombArray(1);
        smartBombStruct.turnX = turnX;
        smartBombStruct.turnY = turnY;
        smartBombStruct.turnCenter = turnCenter;
        smartBombStruct.turnRate = turnRate;
        if ( isempty(smartBombStructs) )
            smartBombStructs = smartBombStruct;
        else
            smartBombStructs(end+1) = smartBombStruct;
        end

    elseif ( isnan(time) )

        playImage = length(smartBombStructs);

    else
        
        goBoomArray = zeros(length(smartBombStructs),2);
        for idxBomb = 1:length(smartBombStructs)
            
            if ( ~isempty(smartBombStructs(idxBomb).pixelsToDelete) )
                playImage = EraseColor(playImage, ...
                    smartBombStructs(idxBomb).pixelsToDelete(1,:),smartBombStructs(idxBomb).pixelsToDelete(2,:),background);
                smartBombStructs(idxBomb).pixelsToDelete = [];
            end
            
            posYx = smartBombStructs(idxBomb).posVector(1:2);
            goBoom = NaN(1,2);
            
            if ( posYx(1) > smartBombStructs(idxBomb).tgtVector(1) )
                goBoom = posYx;
            end
            
            if ( ~isempty(explodeArray) )
                distances = sqrt(sum((explodeArray(2:3,:)-posYx).^2,1));
                if (any((distances-explodeArray(4,:))<=0))
                    goBoom = posYx;
                end
            end
            
            goBoomArray(idxBomb,:) = goBoom;
            
            if ( isnan(goBoom(1)) )
                
                deltaTime = time - smartBombStructs(idxBomb).timeBomb;
                turnAngle = smartBombStructs(idxBomb).turnRate * deltaTime;
                newPos = smartBombStructs(idxBomb).turnCenter + ...
                    smartBombStructs(idxBomb).turnX*cos(turnAngle) + smartBombStructs(idxBomb).turnY*sin(turnAngle);
                newVel = (smartBombStructs(idxBomb).turnY/norm(smartBombStructs(idxBomb).turnY)*cos(turnAngle) - ...
                    smartBombStructs(idxBomb).turnX/norm(smartBombStructs(idxBomb).turnX)*sin(turnAngle));
                destVector = smartBombStructs(idxBomb).tgtVector - newPos;
                explodeRadius = globals.yPixels*globals.explodeRadius;
                
                colours = [0 0 255];
                if ( ~isempty(explodeArray) )
                    
                    relPoss = explodeArray(2:3,:) - newPos(1:2);
                    destDir = destVector(1:2)/norm(destVector(1:2));
                    destDotRel = sum(destDir.*relPoss,1);
                    destDotRel(destDotRel<0) = nan;
                    destDotRel = min(destDotRel,2*explodeRadius);
                    cpaVecs = explodeArray(2:3,:) - (newPos(1:2)+destDir.*destDotRel);
                    cpas = sqrt(cpaVecs(1,:).^2 + cpaVecs(2,:).^2);
                    deltaCpas = cpas - explodeRadius;
                    idx = deltaCpas < 0;
                    if ( any(idx) )
                        destVector = [newPos(1:2) - mean(explodeArray(2:3,idx),2);0];
                        colours = [255 255 0];
                    end
                    
                end
                
                turnAxis = cross(newVel,destVector);
                turnVector = cross(turnAxis,newVel);
                turnRadius = globals.yPixels*globals.smartTurnRadius;
                turnVector = turnRadius*turnVector/norm(turnVector);
                turnCenter = newPos + turnVector;
                turnVectorX = newPos - turnCenter;
                turnVectorY = newVel;
                turnVectorY = norm(turnVectorX)*turnVectorY/norm(turnVectorY);
                
                smartBombStructs(idxBomb).posVector = newPos;
                smartBombStructs(idxBomb).timeBomb = time;
                smartBombStructs(idxBomb).turnX = turnVectorX;
                smartBombStructs(idxBomb).turnY = turnVectorY;
                smartBombStructs(idxBomb).turnCenter = turnCenter;

                yPixels = round(smartBombStructs(idxBomb).posVector(1)) + pixelsToAdd(1,:);
                xPixels = round(smartBombStructs(idxBomb).posVector(2)) + pixelsToAdd(2,:);
                idx = yPixels>=1&yPixels<=globals.yPixels& ...
                    xPixels>=1&xPixels<=globals.xPixels;
                if ( ~isempty(idx) )
                    playImage = SubColor(playImage,yPixels(idx),xPixels(idx),colours);
                    smartBombStructs(idxBomb).pixelsToDelete = [yPixels(idx);xPixels(idx)];
                end
            end
            
        end
        
        idxBoom = find(~isnan(goBoomArray(:,1)))';
        for idxBomb = idxBoom
            if ( ~isnan(goBoomArray(idxBomb,1)) )
                explodeArray(:,end+1) = [time goBoomArray(idxBomb,1) goBoomArray(idxBomb,2) 0];
            end
        end
        smartBombStructs(idxBoom) = [];

    end

end

function [playImage,explodeArray] = HandleCities(playImage,explodeArray,time,background,globals)

    persistent cityStructs cityPixles;

    if ( isinf(time) )
        
        targetLocations = [
            1 + globals.locations(1,:)*(globals.yPixels-1);
            1 + globals.locations(2,:)*(globals.xPixels-1)];
        cityPixles = [
            0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0
            0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0
            0 0 1 1 0 0 0 0 0 0 1 1 0 0 0 0 0 0 0 0 0
            0 0 1 1 0 0 0 1 1 0 1 1 0 0 1 1 0 0 0 0 0
            0 0 1 1 0 0 0 1 1 0 1 1 0 0 1 1 0 0 0 0 0
            0 0 1 1 0 0 0 1 1 0 1 1 0 0 1 1 0 0 0 0 0
            0 0 1 1 0 0 0 1 1 0 1 1 0 0 1 1 0 0 0 0 0
            0 0 1 1 0 0 1 1 1 0 1 1 1 0 1 1 0 0 0 0 0
            0 0 1 1 1 0 1 1 1 0 1 1 1 0 1 1 0 1 1 0 0
            0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 1 1 0 0
            1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1];
        
        yCityCenter = size(cityPixles,1);
        xCityCenter = round(size(cityPixles,2)/2);
        [idx1,idx2] = find(cityPixles);
        cityIndices = [2 3 4 6 7 8];
        for idxCity = 1:length(cityIndices)
            tgtIndex = cityIndices(idxCity);
            posYx = [targetLocations(1,tgtIndex) targetLocations(2,tgtIndex)]';
            yPixels = round(posYx(1)) + idx1 - yCityCenter;
            xPixels = round(posYx(2)) + idx2 - xCityCenter;
            idx = yPixels>=1&yPixels<=globals.yPixels& ...
                xPixels>=1&xPixels<=globals.xPixels;
            if ( ~isempty(idx) )
                playImage = SubColor(playImage,yPixels(idx),xPixels(idx),[255 255 255]);
                cityStruct.cityPixels = [yPixels(idx)';xPixels(idx)'];
            end
            cityStruct.alive = true;
            cityStruct.posYx = posYx;
            if ( isempty(cityStructs) )
                cityStructs = cityStruct;
            else
                cityStructs(idxCity) = cityStruct;
            end
        end
        
    elseif ( time<0 )
        
        for idxCity = 1:length(cityStructs)
            
            playImage = SubColor(playImage, ...
                cityStructs(idxCity).cityPixels(1,:), ...
                cityStructs(idxCity).cityPixels(2,:),[255 255 255]);
            
        end
        
    elseif ( isnan(time) )

        playImage = length(cityStructs);
        
    else

        goBoomArray = zeros(length(cityStructs),2);
        for idxCity = 1:length(cityStructs)

            posYx = cityStructs(idxCity).posYx;
            goBoom = NaN(1,2);

            if ( ~isempty(explodeArray) )
                distances = sqrt(sum((explodeArray(2:3,:)-posYx).^2,1));
                if (any((distances-explodeArray(4,:))<=0))
                    goBoom = posYx;
                end
            end
            
            goBoomArray(idxCity,:) = goBoom;
            
        end
        
        idxBoom = find(~isnan(goBoomArray(:,1)))';
        for idxCity = idxBoom
            if ( ~isnan(goBoomArray(idxCity,1)) )
                explodeArray(:,end+1) = [time goBoomArray(idxCity,1) goBoomArray(idxCity,2) 0];
            end
        end
        cityStructs(idxBoom) = [];

    end
    
end

function [playImage,explodeArray] = HandleExplosions(playImage,explodeArray,time,background,globals)

    persistent explosionTime explosionRadius explosions;
    
    if (time<0)
        
        explosionTime = globals.explodeRadius / globals.explodeSpeed;
        explosionRadius = globals.yPixels * globals.explodeRadius;
        explosions = 0;
        
    else
        
        if (size(explodeArray,2)>explosions)
            soundsc(globals.yExplosion,globals.fExplosion+1000*randn(1));
        end
        explosions = size(explodeArray,2);
        
        yPixels = globals.yPixels;
        xPixels = globals.xPixels;
        removeExplode = zeros(1,size(explodeArray,2))>1;
        for idxExplode = 1:size(explodeArray,2)
            timeSinceExplosion = time - explodeArray(1,idxExplode);
            pointYx = explodeArray(2:3,idxExplode);
            radiusOld = explodeArray(4,idxExplode);
            radiusNew = radiusOld;
            if (timeSinceExplosion<=explosionTime)
                radiusNew = explosionRadius * timeSinceExplosion/explosionTime;
                pixYx = PointRadiusToPixels(yPixels,xPixels,pointYx,radiusNew,radiusOld);
                if (~isempty(pixYx))
                    playImage = AddColor(playImage,pixYx(:,1),pixYx(:,2),[127 127 63]);
                end
            elseif (timeSinceExplosion<=2*explosionTime)
                radiusNew = explosionRadius*(2 - timeSinceExplosion/explosionTime);
                pixYx = PointRadiusToPixels(yPixels,xPixels,pointYx,radiusOld,radiusNew);
                if (~isempty(pixYx))
                    playImage = EraseColor(playImage,pixYx(:,1),pixYx(:,2),background);
                end
            else
                removeExplode(idxExplode) = true;
                pixYx = PointRadiusToPixels(yPixels,xPixels,pointYx,10,0);
                if (~isempty(pixYx))
                    playImage = EraseColor(playImage,pixYx(:,1),pixYx(:,2),background);
                end
            end
            explodeArray(4,idxExplode) = radiusNew;
        end
        explodeArray(:,removeExplode) = [];
        
    end

end

function [playImage,explodeArray] = HandleInterceptor(playImage,explodeArray,time,background,globals)

    persistent interceptorStructs launcherLocations interceptorSpeed ...
        yPixels xPixels lastTime yPixAdd xPixAdd yPixInterceptor xPixInterceptor;

    if (time<0)
        
        targetLocations = [
            1 + globals.locations(1,:)*(globals.yPixels-1);
            1 + globals.locations(2,:)*(globals.xPixels-1)];
        yxSilo = [ 0 0 0 10 10 10 20 20 20;
                  -5 0 5 -5  0  5 -5  0  5];
        yxSilo = [ 0 0 0 10 10 10;
                  -5 0 5 -5  0  5];
        siteLocations = targetLocations(:,[1 5 9]);
        launcherLocations = [];
        for idx = 1:size(yxSilo,2)
            launcherLocations(:,end+1) = siteLocations(:,1) + yxSilo(:,idx);
            launcherLocations(:,end+1) = siteLocations(:,2) + yxSilo(:,idx);
            launcherLocations(:,end+1) = siteLocations(:,3) + yxSilo(:,idx);
        end
        interceptorSpeed = globals.yPixels * globals.interceptorSpeed;
        yPixels = globals.yPixels;
        xPixels = globals.xPixels;
        lastTime = 0;
        yPixAdd = [-3 -2 -1 0  0  0  0 0 0 0 1 2 3]';
        xPixAdd = [ 0  0  0 0 -3 -2 -1 1 2 3 0 0 0]';
        interceptorStructs = [];
        yPixInterceptor = [-3 -2 -1 -1 -1  0 0 0  1 1 1  2 2 2  3 3 3];
        xPixInterceptor = [ 0  0 -1  0  1 -1 0 1 -1 0 1 -1 0 1 -1 0 1];
        for idx = 1:size(launcherLocations,2)
            yPix = round(launcherLocations(1,idx)) + yPixInterceptor;
            xPix = round(launcherLocations(2,idx)) + xPixInterceptor;
            xyIdx = yPix>=1&yPix<=globals.yPixels& ...
                xPix>=1&xPix<=globals.xPixels;
            playImage = SubColor(playImage,yPix(xyIdx),xPix(xyIdx),[255 255 255]);
        end

    elseif (~isempty(globals.controlAxis.UserData))
        
        if (~isempty(launcherLocations))
            
            cursorLocation = [globals.controlAxis.UserData(1,2) globals.controlAxis.UserData(1,1)]';
            yPix = round(cursorLocation(1)) + yPixAdd;
            xPix = round(cursorLocation(2)) + xPixAdd;
            xyIdx = yPix>=1&yPix<=globals.yPixels& ...
                xPix>=1&xPix<=globals.xPixels;
            playImage = SubColor(playImage,yPix(xyIdx),xPix(xyIdx),[255 0 0]);
            
            [~,idx] = min(sum((launcherLocations(2,:)-cursorLocation(2)).^2,1));
            
            startTime = time;
            startYx = launcherLocations(:,idx);
            stopYx = cursorLocation;
            distance = sqrt(sum((stopYx-startYx).^2));
            stopTime = startTime + distance/interceptorSpeed;
            interceptorStruct.yInterpolator = ...
                griddedInterpolant([startTime stopTime],[startYx(1) stopYx(1)]);
            interceptorStruct.xInterpolator = ...
                griddedInterpolant([startTime stopTime],[startYx(2) stopYx(2)]);
            [pixelYx,pixelTimes] = LineToPixels(yPixels,xPixels,startTime,stopTime,startYx,stopYx);
            interceptorStruct.pixelYx = pixelYx;
            interceptorStruct.pixelTimes = pixelTimes;
            if ( isempty(interceptorStructs) )
                interceptorStructs = interceptorStruct;
            else
                interceptorStructs(end+1) = interceptorStruct;
            end
            
            yPix = round(launcherLocations(1,idx)) + yPixInterceptor;
            xPix = round(launcherLocations(2,idx)) + xPixInterceptor;
            xyIdx = yPix>=1&yPix<=globals.yPixels& ...
                xPix>=1&xPix<=globals.xPixels;
            playImage = EraseColor(playImage,yPix(xyIdx),xPix(xyIdx),background);
            launcherLocations(:,idx) = [];
            
            soundsc(globals.yMissile,globals.fMissile+1000*randn(1));
            
        end
        globals.controlAxis.UserData = [];
        
    else

        goBoomArray = zeros(length(interceptorStructs),2);
        for idxInterceptor = 1:length(interceptorStructs)
            
            posYx = [interceptorStructs(idxInterceptor).yInterpolator(time) interceptorStructs(idxInterceptor).xInterpolator(time)]';
            goBoom = NaN(1,2);
            
            if ( time>=interceptorStructs(idxInterceptor).yInterpolator.GridVectors{1}(end) )
                goBoom = [interceptorStructs(idxInterceptor).yInterpolator.Values(end) interceptorStructs(idxInterceptor).xInterpolator.Values(end)];
            end
            
            if ( ~isempty(explodeArray) )
                distances = sqrt(sum((explodeArray(2:3,:)-posYx).^2,1));
                if (any((distances-explodeArray(4,:))<=0))
                    goBoom = posYx;
                end
            end
            
            goBoomArray(idxInterceptor,:) = goBoom;
            
            pixelYx = interceptorStructs(idxInterceptor).pixelYx;
            pixelTimes = interceptorStructs(idxInterceptor).pixelTimes;
            if ( isnan(goBoom(1)) )
                
                idx = pixelTimes>lastTime & pixelTimes<=time;                
                if ( any(idx) )
                    playImage = AddColor(playImage,pixelYx(idx,1),pixelYx(idx,2),[63 63 63]);
                end

            else
                
                idx = pixelTimes<=time;                
                if ( any(idx) )
                    playImage = EraseColor(playImage,pixelYx(idx,1),pixelYx(idx,2),background);
                end
                
            end
            
        end

        idxBoom = find(~isnan(goBoomArray(:,1)))';
        for idxInterceptor = idxBoom
            if ( ~isnan(goBoomArray(idxInterceptor,1)) )
                explodeArray(:,end+1) = [time goBoomArray(idxInterceptor,1) goBoomArray(idxInterceptor,2) 0];
            end
        end
        interceptorStructs(idxBoom) = [];
        
        lastTime = time;
        
        goBoomArray = zeros(size(launcherLocations'));
        for idxLauncher = 1:size(launcherLocations,2)
            
            posYx = launcherLocations(:,idxLauncher);
            goBoom = NaN(1,2);
            
            if ( ~isempty(explodeArray) )
                distances = sqrt(sum((explodeArray(2:3,:)-posYx).^2,1));
                if (any((distances-explodeArray(4,:))<=0))
                    goBoom = posYx;
                end
            end
            goBoomArray(idxLauncher,:) = goBoom;
            
        end

        idxBoom = find(~isnan(goBoomArray(:,1)))';
        for idxLauncher = idxBoom
            if ( ~isnan(goBoomArray(idxLauncher,1)) )
                explodeArray(:,end+1) = [time goBoomArray(idxLauncher,1) goBoomArray(idxLauncher,2) 0];
            end
        end
        launcherLocations(:,idxBoom) = [];
        
    end

end

function [pixelYx,pixelTimes] = LineToPixels(yPixels,xPixels,startTime,stopTime,startYx,stopYx)

    persistent yPix xPix;
    
    if (isempty(yPix))
        yPix = 1:yPixels;
        xPix = 1:xPixels;
    end

    m = (stopYx(1) - startYx(1)) / (stopYx(2) - startYx(2));
    b = startYx(1) - m*startYx(2);
    idx = yPix>=min(startYx(1),stopYx(1)) & yPix<=max(startYx(1),stopYx(1));
    yPixL = yPix(idx);
    idx = xPix>=min(startYx(2),stopYx(2)) & xPix<=max(startYx(2),stopYx(2));
    xPixL = xPix(idx);

    yPixGivenX = m*xPixL + b;
    xPixGivenY = (yPixL-b)/m;
    xyPix = [
        yPixGivenX' xPixL';
        yPixL' xPixGivenY'];
    xyPixFloat = unique(xyPix,'rows');
    xyPixInt = round(xyPixFloat);
    idx = xyPixInt(:,1)>=1 & xyPixInt(:,1)<=yPixels & ...
        xyPixInt(:,2)>=1 & xyPixInt(:,2)<=xPixels;
    xyPixFloat = xyPixFloat(idx,:);
    xyPixInt = xyPixInt(idx,:);
    distance = sqrt(sum((xyPixFloat - startYx').^2,2));
    [distance,idx] = sort(distance);
    xyPixInt = xyPixInt(idx,:);
    totalDistance = sqrt(sum((stopYx-startYx).^2));
    pixTimes = startTime + (stopTime-startTime)*distance/totalDistance;
    idx = pixTimes>=startTime & pixTimes<=stopTime;
    pixelTimes = pixTimes(idx);
    pixelYx = xyPixInt(idx,:);
    
end

function pixelYx = PointRadiusToPixels(yPixels,xPixels,pointYx,radiusB,radiusS)

    persistent yPix xPix;
    
    if (isempty(yPix))
        yPix = 1:yPixels;
        xPix = 1:xPixels;
    end
    
    pixelYx = [];
    
    idxY = yPix >= (pointYx(1)-radiusB) & yPix <= (pointYx(1)+radiusB);
    idxX = xPix >= (pointYx(2)-radiusB) & xPix <= (pointYx(2)+radiusB);
    if (isempty(yPix) | isempty(xPix))
        return;
    end
    
    [yPixGrid,xPixGrid] = meshgrid(yPix(idxY),xPix(idxX));
    distance = sqrt((yPixGrid-pointYx(1)).^2+(xPixGrid-pointYx(2)).^2);
    idx = distance<=radiusB & distance>=radiusS;
    
    yPixGridColumn = reshape(yPixGrid(idx),prod(size(yPixGrid(idx))),1);
    xPixGridColumn = reshape(xPixGrid(idx),prod(size(xPixGrid(idx))),1);
    pixelYx = [yPixGridColumn xPixGridColumn];
    
end

function pixels = SubColor(pixels,i1,i2,colour)

    iz = zeros(size(i1));
    if (isempty(iz))
        return;
    end
    sizePixels = size(pixels);
    i3 = iz+1;
    idx = i1 + (i2-1)*sizePixels(1) + (i3-1)*sizePixels(1)*sizePixels(2);
    pixels(idx) = colour(1);
    i3 = iz+2;
    idx = i1 + (i2-1)*sizePixels(1) + (i3-1)*sizePixels(1)*sizePixels(2);
    pixels(idx) = colour(2);
    i3 = iz+3;
    idx = i1 + (i2-1)*sizePixels(1) + (i3-1)*sizePixels(1)*sizePixels(2);
    pixels(idx) = colour(3);
    
end

function pixels = AddColor(pixels,i1,i2,colour)

    iz = zeros(size(i1));
    if (isempty(iz))
        return;
    end
    sizePixels = size(pixels);
    i3 = iz+1;
    idx = i1 + (i2-1)*sizePixels(1) + (i3-1)*sizePixels(1)*sizePixels(2);
    pixels(idx) = pixels(idx) + colour(1);
    i3 = iz+2;
    idx = i1 + (i2-1)*sizePixels(1) + (i3-1)*sizePixels(1)*sizePixels(2);
    pixels(idx) = pixels(idx) + colour(2);
    i3 = iz+3;
    idx = i1 + (i2-1)*sizePixels(1) + (i3-1)*sizePixels(1)*sizePixels(2);
    pixels(idx) = pixels(idx) + colour(3);
    
end

function pixels = EraseColor(pixels,i1,i2,background)

    iz = zeros(size(i1));
    if (isempty(iz))
        return;
    end
    sizePixels = size(pixels);
    i3 = iz+1;
    idx = i1 + (i2-1)*sizePixels(1) + (i3-1)*sizePixels(1)*sizePixels(2);
    pixels(idx) = background(idx);
    i3 = iz+2;
    idx = i1 + (i2-1)*sizePixels(1) + (i3-1)*sizePixels(1)*sizePixels(2);
    pixels(idx) = background(idx);
    i3 = iz+3;
    idx = i1 + (i2-1)*sizePixels(1) + (i3-1)*sizePixels(1)*sizePixels(2);
    pixels(idx) = background(idx);
    
end

function AxisCallback(src,~)

   src.UserData = get(src, 'CurrentPoint');
   
end

function noise = LaunchNoise()

    dt = 1/8192;
    T = 1;
    t= 0:dt:T;
    f0 = 4000;
    f1 = 2000;
    ydNoise0 = zeros(size(t));
    flipFlop = true;
    timeFlipFlop = 0;
    magFlipFlop = randn(1);
    tend = 1;
    for idx = 1:length(t)
        time = t(idx);
        if (time<=tend)
            freq = f0 + (f1-f0)*time/tend;
            dtf = 1/freq;
            if (time-timeFlipFlop) >= dtf
                flipFlop = ~flipFlop;
                timeFlipFlop = timeFlipFlop+dtf;
                magFlipFlop = randn(1);
            end
            ydNoise0(idx) = magFlipFlop*double(flipFlop);
        end
    end
    reverbTime = 0.07;
    idxReverb = round(reverbTime/dt);
    for kernelidx = 1:5
        reverbKernel(idxReverb*(kernelidx-1)+1) = 1/kernelidx;
    end
    ydNoise = conv(ydNoise0,reverbKernel);
    hissTime = 0.05;
    hissTimeEnd = 0.15;
    ydHiss = zeros(1,round(hissTimeEnd/dt));
    idx = 1:round(hissTime/dt);
    ydHiss(idx) = rand(1,length(idx));
    ydNoise = [ydHiss ydNoise];
    noise = ydNoise;

end

function noise = ExplosionNoise()

    dt = 1/8192;
    T = 2;
    t= 0:dt:T;
    f0 = 1000;
    f1 = 500;
    tend = 1;
    interpPoints = [0 randn(1)];
    time = 0;
    while (time<=tend)
        freq = f0 + (f1-f0)*time/tend;
        dtf = 1/freq;
        time = time+dtf;
        interpPoints(end+1,:) = [time randn(1)];
    end
    interpPoints(end+1,:) = [time+dt 0];
    interpPoints(end+1,:) = [T 0];
    ydNoise0 = interp1(interpPoints(:,1),interpPoints(:,2),t);
    reverbTime = 0.07;
    idxReverb = round(reverbTime/dt);
    for kernelidx = 1:7
        reverbKernel(idxReverb*(kernelidx-1)+1) = 1/kernelidx;
    end
    ydNoise = conv(ydNoise0,reverbKernel);
    noise = ydNoise;

end

function noise = SirenNoise()

    dt = 1/8192;
    T = 2;
    t= 0:dt:T;
    
    for idx = 1:6
        f0 = 200+10*(idx-1)^2;
        f1 = f0+100;
        c = (f1-f0)/T;
        t= 0:dt:T;
        fChirp = 2*pi*(f0*t + c/2*t.^2);
        yChirp(idx,:) = sin(fChirp);
        yHold(idx,:) = sin(2*pi*(f1*t));
    end
    
    sumChirp = sum(yChirp,1);
    sumHold = sum(yHold,1);
    noise = [sumChirp sumHold fliplr(sumChirp)];

end

function [threatSets,colorSchemes] = MakeWaves()

    % 'missile bomber satellite'
    threatSets{1} = { ...
        'time 0 missile 1 speed .05', ...
        'time 0 missile 1 speed .05', ...
        'time 0 missile 1 speed .05', ...
        'time 0 missile 1 speed .05', ...
        'time 5 missile 1 speed .05', ...
        'time 5 missile 1 speed .05', ...
        'time 5 missile 1 speed .05'};
    threatSets{2} = { ...
        'time 0 missile 5 speed .06', ...
        'time 0 missile 3 speed .06', ...
        'time 0 missile 3 speed .06', ...
        'time 0 missile 3 speed .06', ...
        'time 5 missile 5 speed .06', ...
        'time 5 missile 3 speed .06', ...
        'time 5 missile 3 speed .06'};
    threatSets{3} = { ...
        'time 0 missile 5 speed .07', ...
        'time 0 missile 2 speed .07', ...
        'time 0 missile 2 speed .07', ...
        'time 0 missile 2 speed .07', ...
        'time 5 missile 5 speed .07', ...
        'time 5 missile 2 speed .07', ...
        'time 5 missile 2 speed .07', ...
        'time 6 missile 5 speed .07', ...
        'time 7 missile 5 speed .07', ...
        'time 8 missile 5 speed .07', ...
        'time 5 bomber 3 speed .05', ....
        'time 5 bomber 3 speed .05'};
    threatSets{4} = { ...
        'time 0 missile 5 speed .08', ...
        'time 0 missile 2 speed .08', ...
        'time 0 missile 2 speed .08', ...
        'time 0 missile 2 speed .08', ...
        'time 5 missile 5 speed .08', ...
        'time 5 missile 2 speed .08', ...
        'time 5 missile 2 speed .08', ...
        'time 6 missile 3 speed .08', ...
        'time 7 missile 3 speed .08', ...
        'time 8 missile 3 speed .08', ...
        'time 0 bomber 5 speed .06', ...
        'time 0 bomber 3 speed .06', ...
        'time 0 bomber 3 speed .06', ...
        'time 0 bomber 3 speed .06', ...
        'time 5 bomber 5 speed .06', ...
        'time 5 bomber 3 speed .06', ...
        'time 5 bomber 3 speed .06', ...
        'time 5 bomber 3 speed .06'};
     threatSets{5} = { ...
        'time 0 missile 5 speed .09', ...
        'time 0 missile 3 speed .09', ...
        'time 0 missile 3 speed .09', ...
        'time 0 missile 3 speed .09', ...
        'time 5 missile 5 speed .09', ...
        'time 5 missile 3 speed .09', ...
        'time 5 missile 3 speed .09', ...
        'time 6 missile 5 speed .09', ...
        'time 7 missile 5 speed .09', ...
        'time 8 missile 5 speed .09', ...
        'time 0 bomber 5 speed .07', ...
        'time 0 bomber 5 speed .07', ...
        'time 0 bomber 5 speed .07', ...
        'time 2.5 bomber 5 speed .07', ...
        'time 2.5 bomber 5 speed .07', ...
        'time 2.5 bomber 5 speed .07', ...
        'time 5 bomber 5 speed .07', ...
        'time 5 bomber 5 speed .07', ...
        'time 5 bomber 5 speed .07', ...
        'time 7.5 bomber 5 speed .07', ...
        'time 7.5 bomber 5 speed .07', ...
        'time 7.5 bomber 5 speed .07', ...
        'time 10 bomber 5 speed .07' ...
        'time 10 bomber 5 speed .07' ...
        'time 10 bomber 5 speed .07'};
    
    colorSchemes{1} = [127 0 0;
                       63 0 0];
    colorSchemes{2} = [255 127 0;
                       127 63 0];
    colorSchemes{3} = [127 191 255;
                       32 223 63];
    colorSchemes{4} = [255 63 0;
                       127 0 0];
    colorSchemes{5} = [0 0 63;
                       0 0 31];                   
                   
end
