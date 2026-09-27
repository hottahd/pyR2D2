"""
Functions for push information to Google Spreadsheet
"""

import glob
import os


def init_gspread(json_key, project):
    """
    This function initialize the utility of google spread

    Parameters
    ----------
    json_key : str
        file of json key to access Google API
    project : str
        Project name, typically name of upper directory

    Returns
    -------
    gc : gspread.client.Client
        Instance of Google API
    """

    import gspread
    from google.oauth2.service_account import Credentials

    scopes = [
        "https://spreadsheets.google.com/feeds",
        "https://www.googleapis.com/auth/drive",
    ]

    credentials = Credentials.from_service_account_file(json_key, scopes=scopes)
    gc = gspread.authorize(credentials)

    return gc


################################################################################
def fetch_URL_gspread(json_key=None, project=None):
    """
    Fetchs corresponding URL for google spreadsheet

    Parameters
    ----------
    json_key : str
        File of json key to access Google API
    projet : str
        Project name, typically name of upper directory

    Returns
    -------
    URL : str
        Google spreadsheet URL
    """
    if project == None:
        project = os.getcwd().split("/")[-2]

    if json_key == None:
        json_key = glob.glob(os.environ["HOME"] + "/json/*")[0]
    gid = init_gspread(json_key, project).open(project).id

    return "https://docs.google.com/spreadsheets/d/" + gid


################################################################################


def set_top_line(json_key=None, project=None):
    """
    This function set top line of google spreadsheet

    台帳の管理列 (:py:mod:`pyR2D2.ledger`) のうちシートに無い見出しを、
    既存の見出しの右端に足す。既存の列 (人が書く列を含む) は動かさない。

    Parameters
    ----------
    json_key : str
        File of json key to access Google API
    projet : str
        Project name, typically name of upper directory

    Returns
    -------
        None
    """
    from pyR2D2.ledger import collect, sheet

    if project == None:
        project = os.getcwd().split("/")[-2]

    if json_key == None:
        json_key = glob.glob(os.environ["HOME"] + "/json/*")[0]
    gc = init_gspread(json_key, project)
    wks = gc.open(project).sheet1

    plan = sheet.plan_updates(wks.get_all_values(), [], collect.default_server())
    sheet.apply_plan(wks, plan)


################################################################################
def set_cells_gspread(data, json_key=None, project=None, caseid=None):
    """
    Outputs parameters to Google spreadsheet

    ``r2d2plus-ledger push`` と同じ処理 (:py:func:`pyR2D2.ledger.cli.push_records`)
    で、そのランの管理列のセルだけを見出しの名前で探して書く。人が書く列には
    触らない (以前は A〜T 列を位置で上書きし、21 個目の Server が落ちていた)。

    Parameters
    ----------
    data : pyR2D2.Data, or, pyR2D2.Read
        instance of pyR2D2.Data or pyR2D2.Read classes
    json_key : str
        File of json key to access Google API
    project : str
        Project name, typically name of upper directory
    caseid : str
        Case ID

    """
    from pathlib import Path

    from pyR2D2.ledger import cli, collect

    if project == None:
        project = os.getcwd().split("/")[-2]

    if json_key == None:
        json_key = glob.glob(os.environ["HOME"] + "/json/*")[0]

    run_dir = Path(data.datadir).resolve().parent
    server = collect.default_server()
    record = collect.collect_run(run_dir, server=server)
    if caseid is not None:
        record.caseid = caseid
        record.values["Case ID"] = caseid

    gc = init_gspread(json_key, project)
    wks = gc.open(project).sheet1
    cli.push_records([record], wks, server, project)
