
class SWSimulation {

    // class constructor
    /**
     * Sample JSDoc comment
     */
    constructor() {
        console.log("Initializing SWSimulation class.");

        var self = this;
        self.is_running = false;
        self.n_steps_run = 0;

        // Menu of field types with their label codes
        self.fields = [
            { val: 0, text: "UUGF" },
            { val: 1, text: "VVGF" }
            //{ val: 0, text: "UUGF" }

        ];
        // Default field to plot '0'
        self.fieldPlotted = 0;

        $('#control #play').on('click', function () {
            console.log("Setting simulation to running.");
            self.is_running = true;
        });

        $('#control #pause').on('click', function () {
            console.log("Setting simulation to paused.");
            self.is_running = false;
        });

        var select = $('<select>');
        $(self.fields).each(function () {
            select.append($("<option>").attr('value', this.val).text(this.text));
        });
        select.change(function () {
            var sel = this;
            self.fieldPlotted = sel.value;
            self.updatePlotData(sel.value);
        });
        $('#field').append(select);

        Module._init_sim();
        Module._run_sim(1);
        self.showPlot(self.fieldPlotted);
        self.runSim();
    }

    runSim() {
        var self = this;
        async_recurse_with_test(100, self.runStep, { 'self': self }, self.isRunning, { 'self': self });
    }

    runStep(_self) {
        // Iterate next step of the simulation
        var self = _self['self'];

        Module._run_sim(1);
        // Add one to counter
        self.n_steps_run++;
        $('#counter').text(self.n_steps_run);
        self.updatePlotData(self.fieldPlotted);
    }

    isRunning(_self) {
        return _self['self'].is_running;
    }

    showPlot(g) {
        // This sets the plot type, settings and data input.
        var self = this;
        var arr = get_SW_arrays(g);

        // Rework for ScalarWave GFS
        var data = [{
            //colorscale : [['0' , 'rgb(0,225,100)'],['1', 'rgb(100,0,200)']]
            x: arr['x'],
            y: arr['y'],
            a: arr['a'],
            b: arr['b'],
            type: 'carpet',
            aaxis: {
                smoothing: 0,
                minorgridcount: 9,
                type: 'linear',
                showticklabels: "none",
                showgrid: false // Original was false; wanted to see if true sheds light on graphing issue
            },
            baxis: { // radially outward axis
                smoothing: 0,
                minorgridcount: 9,
                type: 'linear',
                showticklabels: "none",
                showgrid: false // Original was false; wanted to see if true sheds light on graphing issue
            }
        }, {
            z: arr['z'],
            a: arr['a'],
            b: arr['b'],
            type: 'contourcarpet'
        }];
        // Needs to be  3-D
        var layout = {
            yaxis: {
                zeroline: true,
                type: 'linear',
                showgrid: true, // Original was false
                autorange: true,
                showline: false,
                ticks: '',
                showticklabels: true // Original was false, gives numerical labels to y-axis grid
            },
            xaxis: {
                zeroline: true,
                type: 'linear',
                scaleratio: 1,
                scaleanchor: "y", // Scales x-axis to the size of the y-axis
                showgrid: true, // Original was false
                autorange: true,
                showline: false,
                ticks: '',
                showticklabels: true // Original was false, gives numerical labels to x-axis grid
            },
            zaxis: {
                zeroline: true,
                type: 'linear',
                scaleratio: 1,
                scaleanchor: "y", // Scales z-axis to the size of the y-axis
                showgrid: true, // Original was false
                autorange: true,
                showline: false,
                ticks: '',
                showticklabels: true, // Original was false, gives numerical labels to x-axis grid
            },
            hovermode: "closest",
            height: 900 // Original was 700 and was slightly smaller than desired
        };


        Plotly.newPlot('display', data, layout);
        Plotly.Plots.resize('display');
        //Plotly.newPlot('myDiv', data, layout, {displaylogo: false});
    }

    updatePlotData(g) {
        console.log('Updating plot data, showing field', g, '.');
        var arr = get_SW_arrays(g);
        // vv this will probably be redundant 
        //Plotly.restyle('display', { z: [arr['z']] });
    }

} // SWSimulation class



function sleep(ms) {
    return new Promise(resolve => setTimeout(resolve, ms));
}

async function async_recurse_with_test(sleep_time, call, call_args, test, test_args) {
    if (test(test_args)) {
        call(call_args);
    }

    await sleep(sleep_time);
    async_recurse_with_test(sleep_time, call, call_args, test, test_args);
}

// Leave this for now till variables are understood
function get_SW_arrays(g) {

    // TODO: redo all this neatly.
    // 
    var gfs_ptr_F64 = Module._get_gfs() / 9; // /8 does this need to be divided? >> two grid functions incoming
    var gfs_size = Module._get_gfs_size();
    var gfs_arr = Module.HEAPF64.subarray(gfs_ptr_F64, gfs_ptr_F64 + gfs_size);


    
    // split incoming gfs into 
    var nx = 128;
    var ny = 32;
    var nz = 2;

    var NGHOSTS = 3;
    var Npts1 = nx + 2 * NGHOSTS;
    var Npts2 = ny + 2 * NGHOSTS;
    var Npts3 = nz + 2 * NGHOSTS;

    var x_arr = new Array(2 * nx * ny);
    var y_arr = new Array(2 * nx * ny);
    var z_arr = new Array(2 * nx * ny);
    var a_arr = new Array(2 * nx * ny);
    var b_arr = new Array(2 * nx * ny);

    for (var i = 0; i < nx; i++) {
        for (var j = 0; j < ny; j++) {
            var k = 0;
            var g_idx = (i + NGHOSTS) + Npts1 * ((j + NGHOSTS) + Npts2 * ((k + NGHOSTS) + Npts3 * g));
            var a_idx = i * ny + j;
            var a_idx_refl = nx * ny + i * ny + j;

            var r = (i + 1) / nx; //(i+1.0)
            var theta = (j + 0.5) * 1.0 / ny * Math.PI; //(j+1.0)
            var x_of_pt = r * Math.sin(theta);
            var y_of_pt = r * Math.cos(theta);

            a_arr[a_idx] = i;
            b_arr[a_idx] = j;
            x_arr[a_idx] = x_of_pt;
            y_arr[a_idx] = y_of_pt;
            z_arr[a_idx] = gfs_arr[g_idx];
            // Is this where the plot mirroring goes wrong?
            b_arr[a_idx_refl] = -j; // Default is -j
            a_arr[a_idx_refl] = -i;
            y_arr[a_idx_refl] = -y_of_pt;
            x_arr[a_idx_refl] = -x_of_pt;

            z_arr[a_idx_refl] = gfs_arr[g_idx];
        }
    }

    return {
        a: a_arr,
        b: b_arr,
        x: x_arr,
        y: y_arr,
        z: z_arr,
    };
}

