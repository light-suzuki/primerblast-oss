"""Thread-local, optional workflow observation; no change to PCR decisions."""
from contextlib import contextmanager
from contextvars import ContextVar
import time

_observer = ContextVar('workflow_observer', default=None)


class ProgressRecorder:
    def __init__(self, callback):
        self.callback = callback
        self.started = time.monotonic()
        self.phase_started = self.started
        self.stage = None
        self.context = {}
        self.counts = {}
        self.seconds = {}

    def report(self, stage, completed=None, total=None, **context):
        now = time.monotonic()
        if stage != self.stage:
            if self.stage is not None:
                self.seconds[self.stage] = self.seconds.get(self.stage, 0) + now - self.phase_started
            self.phase_started = now
            self.stage = stage
            self.counts = {}
        self.context.update(context)
        if completed is not None:
            self.counts['completed'] = completed
        if total is not None:
            self.counts['total'] = total
        self.callback({'progress': dict(self.context, **self.counts,
            stage=stage, elapsed_seconds=now - self.started,
            stage_elapsed_seconds=now - self.phase_started,
            stage_seconds=dict(self.seconds))})

    def finish(self, stage='done'):
        self.report(stage)
        return {'total_seconds': time.monotonic() - self.started,
                'stage_seconds': dict(self.seconds)}


@contextmanager
def observe(callback):
    recorder = ProgressRecorder(callback)
    token = _observer.set(recorder)
    try:
        yield recorder
    finally:
        _observer.reset(token)


def report(stage, completed=None, total=None, **context):
    recorder = _observer.get()
    if recorder is not None:
        recorder.report(stage, completed, total, **context)


def partial(result):
    recorder = _observer.get()
    if recorder is not None:
        recorder.callback({'partial_result': result})


def active():
    return _observer.get() is not None
