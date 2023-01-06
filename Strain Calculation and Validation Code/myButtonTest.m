function myButtonTest(fig,cbx1,cbx2,cbx3,cbx4,cbx5,cbx6,cbx7,cbxY)
PushButton = uicontrol(fig,'Style', 'push', 'String', 'Ok','Position', [20 20 100 20],'CallBack', @PushB);
uiwait(fig)
    function PushB(source,event)
        global results
        results = [cbx1.Value cbx2.Value cbx3.Value cbx4.Value cbx5.Value cbx6.Value cbx7.Value cbxY.Value];
        close(fig)
    end
end