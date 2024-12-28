from typing import Literal


class TaskManager:
    def __init__(self, tasks = []):
        self._tasks = tasks

    def add_task(self, task):
        self._tasks.append(task)
        return

    @property
    def all_tasks(self):
        return self._tasks

    def get_overview(self):
        return [task.description for task in self._tasks]


STATUS_OPTIONS = Literal['r', 'rb', 'w', 'wb']

class task:
    def __init__(self, description:str, status = 'not done'):
        self.description = description
        self._status = status

    def set_stat8s(self, status: STATUS_OPTIONS):
        pass

    def set_complete(self):
        self._status = 'complete'

    def set_in_progress(self):
        self._status = 'in progress'

    def show_status(self):
        return self._status

task1 = task('do laundry')
task2 = task('parapy interview')

to_do = TaskManager([task1, task2])

print(to_do.get_overview())

task3 = task('foo')

to_do.add_task(task3)

to_do._tasks

task1.set_in_progress()

print(task1.set_stat8s('remco'))


class Remco:
    def __init__(self):
        pass

    colleague = "max"



