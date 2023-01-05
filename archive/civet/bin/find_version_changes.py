#! /usr/bin/env python

"""
The civet_research_pipelines repo contains many pipelines.  In any particular
new version of the released repo, most pipelines are unaffected.  But when
were the pipelines changed?  This program attempts to answer the question for
each pipeline in the repo.

This program depends on a utility written by Al named find_duplicate_files.py.
It creates a database of file information including MD5 checksums for every
file in a directory tree.  Before running this program, find_duplicate_files.py
must be run over the entire /opt/compsci/civet_research_pipelines tree, with the
command line:

find_duplicate_files.py -d civet_research_files.db -m 0 -s .git /opt/compsci/civet_research_pipelines

The basic outline for this program is for every file in the database:
 - split the full file path in three parts:
   - /opt/compsci/civet_research_pipelines
   - version number
   - pipeline-and-file
 - isolate the pipeline from pipeline-and-file, where pipeline is everything
   before the final filename OR bin/filename
 - for each pipeline -> file -> md5sum record the earliest version.
 - for each pipeline build a list of all of its files' earliest versions.  Those
   are the versions in which the pipeline changed.

The schema for the database is:

sqlite> .schema
CREATE TABLE files (
        id INTEGER PRIMARY KEY,
        path TEXT KEY NOT NULL,
        size INTEGER KEY NOT NULL,
        mtime INTEGER NOT NULL,
        inode INTEGER NOT NULL,
        small_md5 TEXT KEY,
        large_md5 TEXT KEY,
        seen_this_pass INTEGER KEY NOT NULL,
        UNIQUE (path, size, mtime, inode));

Files which do not have a duplicated small_md5 will have large_md5 == None. If
that happens the file is unique to that version, so record it.

"""
import argparse
import sqlalchemy
from sqlalchemy import Column, Integer, String
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker

from pipeline_definitions import pipes

def parse_args():
    parser = argparse.ArgumentParser("find version changes V0")
    parser.add_argument(
        'database',
        help="The path to the database of civet_research_pipeline files"
    )
    parser.add_argument('--by-version', action='store_true',
                        help="Return the list of pipelines that changed in a "
                             "version, instead of the list of versions in "
                             "which a pipeline changed.")
    parser.add_argument('--simple-names', action='store_true',
                        help="Use short pipeline command names instead of "
                             "directory paths.")
    return parser.parse_args()

Base = declarative_base()


class CRPFile(Base):
    __tablename__ = 'files'
    id = Column(Integer, primary_key=True)
    path = Column(String(250), nullable=False)
    size = Column(Integer, nullable=False)
    mtime = Column(Integer, nullable=False)
    inode = Column(Integer, nullable=False)
    small_md5 = Column(String(32))
    large_md5 = Column(String(32))
    seen_this_pass = Column(Integer, nullable=False)

    def extract(self):
        parts = self.path.split('/')
        version = parts[4]
        if parts[-2] == 'bin':
            path_parts = parts[5:-2]
            filename = parts[-2:]
        else:
            path_parts = parts[5:-1]
            filename = parts[-1:]

        return version, path_parts, filename, self.large_md5


def cmp_to_key(mycmp):
    """
    This is hacked from functools.cmp_to_key()
    Convert a cmp= function into a key= function'
    The real one has cases for all possible numeric comparisons, but the
    sorted() function only uses less-than.
    We use this as the lambda to sort version numbers.
    :param mycmp: A comparison routine.
    :return: True if self is less than other.
    """
    class K(object):
        def __init__(self, obj, *args):
            self.obj = obj
        def __lt__(self, other):
            return not mycmp(self.obj, other.obj)
    return K


def get_engine(db):
    return sqlalchemy.create_engine('sqlite:///' + db)


pipe_changes = {}


def get_name_dict():
    name_dict = {}
    for cmd in pipes.keys():
        name_dict[pipes[cmd]['dir']] = cmd
    return name_dict


def print_by_version(simple_names):
    """
    Walk the nested dict of pipe_changes, and collect all the pipelines (paths)
    that changed in each version.
    :param simple_names: Use short command names  (if known) instead of
        pipeline paths.
    :return: None
    """
    versions = {}
    if simple_names:
        name_dict = get_name_dict()
    paths = sorted(list(pipe_changes.keys()))
    for p in paths:
        pipe = pipe_changes[p]
        for f in pipe.keys():
            filename = pipe[f]
            for m in filename.keys():
                try:
                    # filename[m] is a version number. p is a pipeline path.
                    versions[filename[m]].add(p)
                except KeyError:
                    versions[filename[m]] = {p}

    for v in sorted(list(versions.keys()),  key=cmp_to_key(new_is_earlier)):
        changed_pipes = versions[v]
        print(v)
        for p in sorted(list(changed_pipes)):
            if simple_names:
                # Look up the path and return the simple command name, if the
                # path is known to have a simple command. Else, return the path.
                p = name_dict.get(p, p)

            print('   ', p)


def print_by_pipeline(simple_names):
    """
    Walk the nested dict of pipe_changes, and collect all the lowest versions
    recorded for any file in a path into a set.
    Each lowest version marks a change in that pipeline.
    :param simple_names: Use short command names  (if known) instead of
        pipeline paths.
    :return: None
    """
    if simple_names:
        name_dict = get_name_dict()
    paths = sorted(list(pipe_changes.keys()))
    for p in paths:
        changes = set()
        pipe = pipe_changes[p]
        for f in pipe.keys():
            filename = pipe[f]
            for m in filename.keys():
                changes.add(filename[m])
        if simple_names:
            # Look up the path and return the simple command name, if the path
            # is known to have a simple command. Else, return the path.
            p = name_dict.get(p, p)

        print(p, ', '.join(sorted(list(changes),
                                  key=cmp_to_key(new_is_earlier))))


def new_is_earlier(cv, nv, path, filename):
    """
    These are dot-separated version numbers.
    Go toward  more minor numbers until there is a difference. If one is longer
    than the other, and all of the fields in the shorter one match, then the
    longer one is later.  Remember to treat all as ints.
    :param cv: the current version string.
    :param nv: the new version string, which may be earlier than the current
        one.
    :param path: the path to the file whose version is being checked.
        Used only for error reporting.
    :param filename: the name of the file whose version is being checked.
        Used only for error reporting.
    :return: True iff the new version string represents an earlier version.
    """
    try:
        cv_parts = [int(x) for x in cv.split('.')]
        nv_parts = [int(x) for x in nv.split('.')]
    except:
        # Catch non-integer version values
        return False
    for n in range(len(cv_parts)):
        try:
            if nv_parts[n] == cv_parts[n]:
                continue
            return nv_parts[n] < cv_parts[n]
        except IndexError:
            # They have matched up to here, but cv_parts is longer than
            # nv_parts. Therefore, nv is earlier.
            return True

    # If we get here, all the elements up to the length of the shorter version
    # string are equal. We have one of two situations. Either they are the same
    # length and completely equal (not sure how that happens, but...), or
    # nv_parts is longer than cv_parts. For both situations, nv is not newer.
    if cv == nv:
        print("Two files with the same version number but different "
              "checksums:\n"
              "    {}/{} ({})".format(path, filename, cv))
    return False


def add_lowest_version(path, filename, md5, version):
    """
    Build a nested dictionary containing the earliest version for which a
    particular md5 of a file was found.  This will result in one entry for
    each version in which a file changed.  When all the change versions for all
    the files in a pipeline are collected, we have the set of versions in which
    the pipeline changed.
    :param path: the path to the file: this is usually the pipeline.
    :param filename: the filename that we are checking.
    :param md5: the md5sum of the file
    :param version: the version containing this file.
    :return: None
    """
    try:
        path_dict = pipe_changes[path]
    except KeyError:
        path_dict = {}
        pipe_changes[path] = path_dict
    try:
        file_dict = path_dict[filename]
    except KeyError:
        file_dict = {}
        path_dict[filename] = file_dict
    try:
        # The file dict contains the minimum version for each
        current_minimum_version = file_dict[md5]
        if new_is_earlier(current_minimum_version, version):
            file_dict[md5] = version
    except KeyError:
        # Hadn't seen this md5 for this file before.  So this is the earliest
        file_dict[md5] = version


def main():
    args = parse_args()
    engine = get_engine(args.database)

    DBSession = sessionmaker(bind=engine)
    session = DBSession()

    crp_files = session.query(CRPFile)
    for crp_file in crp_files:
        version, path, filename, md5 = crp_file.extract()
        if not path:
            continue
        # Only handle numeric versions, not test versions.
        if not version[0].isdigit():
            continue

        add_lowest_version('/'.join(path), '/'.join(filename), md5, version)

    if args.by_version:
        print_by_version(args.simple_names)
    else:
        print_by_pipeline(args.simple_names)

if __name__ == '__main__':
    main()
