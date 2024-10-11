# trace_calls.py

import gdb

class TraceAllCalls(gdb.Command):
    """Trace all function calls."""

    def __init__(self):
        super(TraceAllCalls, self).__init__("tracecalls", gdb.COMMAND_USER)
        self.enabled = False

    def invoke(self, arg, from_tty):
        if not self.enabled:
            self.enabled = True
            gdb.events.stop.connect(self.on_stop)
            print("Function call tracing enabled.")
        else:
            self.enabled = False
            gdb.events.stop.disconnect(self.on_stop)
            print("Function call tracing disabled.")

    def on_stop(self, event):
        frame = gdb.newest_frame()
        if frame:
            func_name = frame.name()
            if func_name:
                print("Function call: {}".format(func_name))
        gdb.execute('continue', to_string=True)

TraceAllCalls()
