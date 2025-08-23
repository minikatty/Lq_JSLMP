function switchOccurred = detectSwitch(currentState)
    persistent prevState
    persistent hasTriggered
    
    % 初始化
    if isempty(prevState)
        prevState = currentState;
        hasTriggered = false;
        switchOccurred = false;
        return;
    end
    
    % 检测A→B切换，且只触发一次
    if strcmp(prevState, 'Adam') && strcmp(currentState, 'VSS') && ~hasTriggered
        switchOccurred = true;
        hasTriggered = true;
        disp('检测到步长机制切换！');
    else
        switchOccurred = false;
    end
    
    % 当状态不再是B时，重置触发标志
    if ~strcmp(currentState, 'VSS')
        hasTriggered = false;
    end
    
    prevState = currentState;
end