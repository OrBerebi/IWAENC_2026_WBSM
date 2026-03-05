% Helper function triggered when the loop is on and the player stops.
% Player is restarted to continue loop.
function player_stop(obj, event)
    stop(obj);
    play(obj);